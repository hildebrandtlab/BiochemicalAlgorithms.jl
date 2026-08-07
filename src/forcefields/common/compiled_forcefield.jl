export
    CompiledForceField,
    compile,
    finalize_optimization!,
    rebuild_pairlist!,
    shuffle_interactions!

# ============================================================================
#  Cached, struct-of-arrays (SoA), index-based evaluator for a `ForceField`.
#
#  Rationale: the reference path in `optimize_structure!` rebuilds
#  the entire nonbonded interaction list on *every* energy evaluation and
#  accesses atom coordinates/forces through `Atom{T}` row-views (each a
#  `Dict{Int,Int}` lookup). This evaluator lowers the already-correct output of
#  `setup!`/`update!` into flat arrays indexed by *row number* into contiguous
#  coordinate/force buffers, so per-evaluation cost is a tight allocation-free
#  loop. The pair list is rebuilt only on demand (BALL-style), not every eval.
#
#  The existing per-component `compute_energy!`/`compute_forces!`/`update!`
#  remain the slow *reference* path and are untouched; this is purely additive.
#
#  All energy/force formulas below are transcribed VERBATIM from
#  stretch_component.jl / bend_component.jl / torsion_component.jl /
#  nonbonded_component.jl so results match the reference path bit-for-bit when
#  the accumulation type `Acc` equals the storage type `T`.
# ============================================================================

# --- evaluation backends -----------------------------------------------------
abstract type EvalBackend end
struct SerialBackend <: EvalBackend end
struct ThreadedBackend <: EvalBackend
    nthreads::Int
end

# --- per-component SoA payloads (row numbers into r/F buffers) ---------------
struct StretchArrays{Acc}
    i::Vector{Int32}
    j::Vector{Int32}
    r0::Vector{Acc}
    k::Vector{Acc}
end
StretchArrays{Acc}() where Acc = StretchArrays{Acc}(Int32[], Int32[], Acc[], Acc[])

struct BendArrays{Acc}
    i::Vector{Int32}   # a1
    j::Vector{Int32}   # a2 (center)
    l::Vector{Int32}   # a3
    θ0::Vector{Acc}
    k::Vector{Acc}
end
BendArrays{Acc}() where Acc = BendArrays{Acc}(Int32[], Int32[], Int32[], Acc[], Acc[])

# Torsions have a variable number of Fourier terms; store CSR-flattened.
struct TorsionArrays{Acc}
    i::Vector{Int32}
    j::Vector{Int32}
    k::Vector{Int32}
    l::Vector{Int32}
    term_offset::Vector{Int32}   # length nt+1, CSR offsets into V/ϕ0/f/div
    V::Vector{Acc}
    ϕ0::Vector{Acc}
    f::Vector{Int32}
    div::Vector{Int32}
end
TorsionArrays{Acc}() where Acc =
    TorsionArrays{Acc}(Int32[], Int32[], Int32[], Int32[], Int32[1], Acc[], Acc[], Int32[], Int32[])

# vdW (12-6) and hydrogen bonds (12-10) share this layout; `es_*` carries the
# matching electrostatic contribution so a pair contributes both in one pass is
# *not* assumed here — electrostatics are kept in their own array to mirror the
# reference exactly.
struct PairLJArrays{Acc}
    i::Vector{Int32}
    j::Vector{Int32}
    A::Vector{Acc}
    B::Vector{Acc}
    scaling::Vector{Acc}
end
PairLJArrays{Acc}() where Acc = PairLJArrays{Acc}(Int32[], Int32[], Acc[], Acc[], Acc[])

struct PairESArrays{Acc}
    i::Vector{Int32}
    j::Vector{Int32}
    q1q2::Vector{Acc}
    scaling::Vector{Acc}
end
PairESArrays{Acc}() where Acc = PairESArrays{Acc}(Int32[], Int32[], Acc[], Acc[])

# --- the compiled evaluator --------------------------------------------------
mutable struct CompiledForceField{T,Acc,B<:EvalBackend}
    ff::ForceField{T}
    backend::B

    # bonded (never rebuilt during minimization)
    stretch::StretchArrays{Acc}
    bend::BendArrays{Acc}
    proper::TorsionArrays{Acc}
    improper::TorsionArrays{Acc}

    # nonbonded (rebuilt on demand)
    lj::PairLJArrays{Acc}
    hb::PairLJArrays{Acc}
    es::PairESArrays{Acc}

    # cached nonbonded scalars
    vdw_sw::CubicSwitchingFunction{Acc}
    es_sw::CubicSwitchingFunction{Acc}
    es_prefactor::Acc
    es_prefactor_force::Acc
    distance_dependent_dielectric::Bool

    # coordinate / force scratch (row-indexed, Acc precision)
    r::Vector{Vector3{Acc}}
    F::Vector{Vector3{Acc}}

    # per-thread scratch for the ThreadedBackend (empty for SerialBackend)
    F_threads::Vector{Vector{Vector3{Acc}}}
    e_threads::Vector{Acc}

    # per-category energies (keys match ForceField reference path)
    energy::Dict{String,Acc}

    # idx <-> row bookkeeping
    idx2row::Dict{Int,Int32}
    natoms::Int

    # constrained atoms, as a per-row mask
    constrained_mask::Vector{Bool}

    # BALL-style pair-list rebuild policy
    update_frequency::Int
    max_displacement::Acc
    r_at_last_build::Vector{Vector3{Acc}}
    iters_since_build::Int

    # device-resident state for the GPU backend (nothing for CPU backends)
    gpu::Any
end

@inline accumulation_type(::CompiledForceField{T,Acc}) where {T,Acc} = Acc

function Base.show(io::IO, cff::CompiledForceField{T,Acc,B}) where {T,Acc,B}
    print(io, "CompiledForceField{$T, acc=$Acc, $(nameof(B))}: ",
        cff.natoms, " atoms, ",
        length(cff.stretch.k), " stretches, ",
        length(cff.bend.k), " bends, ",
        length(cff.proper.i), "+", length(cff.improper.i), " torsions, ",
        length(cff.lj.i), " vdW / ", length(cff.hb.i), " hbond / ",
        length(cff.es.i), " electrostatic pairs")
end

# ----------------------------------------------------------------------------
#  Compilation: lower the reference components into SoA arrays.
# ----------------------------------------------------------------------------
"""
    compile(ff::ForceField{T}; acc=T, backend=:serial, update_frequency=1, max_displacement=8)

Build a [`CompiledForceField`](@ref) from `ff`: a cached, struct-of-arrays
evaluator suitable for repeated energy/force evaluation (e.g. inside
`optimize_structure!`). `acc` selects the accumulation precision (`Float32` or
`Float64`); `backend` selects `:serial` (more backends added in later stages).

Calls `update!(ff)` once to obtain the correct interaction list, then lowers it.
The atom topology (and therefore row numbering) must not change for the lifetime
of the returned object.
"""
function compile(ff::ForceField{T};
    acc::Type{Acc}=T,
    backend::Symbol=:serial,
    update_frequency::Integer=0,
    max_displacement::Real=-1) where {T,Acc<:Real}

    b = _select_backend(Val(backend), ff, Acc)

    at = atoms(ff.system)
    natoms = length(at.idx)
    idx2row = Dict{Int,Int32}(idx => Int32(k) for (k, idx) in enumerate(at.idx))

    r = Vector{Vector3{Acc}}(undef, natoms)
    F = Vector{Vector3{Acc}}(undef, natoms)

    nthreads = b isa ThreadedBackend ? b.nthreads : 0
    F_threads = [Vector{Vector3{Acc}}(undef, natoms) for _ in 1:nthreads]
    e_threads = zeros(Acc, nthreads)

    constrained_mask = falses(natoms)
    for a in ff.constrained_atoms
        # `constrained_atoms` holds row positions into atoms(system); mirror the
        # reference path which indexes atoms(ff.system)[constrained_atoms].
        if 1 <= a <= natoms
            constrained_mask[a] = true
        end
    end

    cff = CompiledForceField{T,Acc,typeof(b)}(
        ff, b,
        StretchArrays{Acc}(), BendArrays{Acc}(), TorsionArrays{Acc}(), TorsionArrays{Acc}(),
        PairLJArrays{Acc}(), PairLJArrays{Acc}(), PairESArrays{Acc}(),
        CubicSwitchingFunction{Acc}(zero(Acc), zero(Acc)),
        CubicSwitchingFunction{Acc}(zero(Acc), zero(Acc)),
        Acc(ES_Prefactor), Acc(ES_Prefactor_force), false,
        r, F, F_threads, e_threads,
        Dict{String,Acc}(),
        idx2row, natoms,
        constrained_mask,
        Int(update_frequency), Acc(max_displacement),
        Vector{Vector3{Acc}}(undef, natoms), 0,
        nothing,
    )

    # auto-derive a safe Verlet skin from the cutoffs unless overridden: a pair
    # absent from the list was >= nonbonded_cutoff apart at build time; it can
    # only reach the (vdw/es) cutoff if atoms move > (nonbonded_cutoff-cutoff)/2.
    if max_displacement < 0
        nb = Acc(ff.options[:nonbonded_cutoff])
        cut = max(Acc(ff.options[:vdw_cutoff]), Acc(ff.options[:electrostatic_cutoff]))
        cff.max_displacement = max(zero(Acc), (nb - cut) / 2)
    end

    _lower_bonded!(cff)
    sync_positions_from_table!(cff)   # fill cff.r from the table FIRST
    rebuild_pairlist!(cff)            # then lower nonbonded + snapshot positions
    _backend_init!(cff)               # GPU extensions allocate/upload device state here
    cff
end

_select_backend(::Val{:serial}, ::ForceField, ::Type) = SerialBackend()
_select_backend(::Val{:threads}, ::ForceField, ::Type) = ThreadedBackend(Threads.nthreads())
function _select_backend(::Val{:gpu}, ff::ForceField, ::Type{Acc}) where Acc
    hasmethod(_make_gpu_backend, Tuple{ForceField,Type}) || throw(ArgumentError(
        "the :gpu backend needs `using KernelAbstractions` plus a device package " *
        "(`using Metal` on Apple Silicon, or `using CUDA`)"))
    _make_gpu_backend(ff, Acc)
end
_select_backend(::Val{B}, ::ForceField, ::Type) where B =
    throw(ArgumentError("unknown backend :$B (available: :serial, :threads, :gpu)"))

# Implemented by the KernelAbstractions extension (BiochemicalAlgorithmsGPUExt).
function _make_gpu_backend end

# GPU device registry. Filled by the Metal/CUDA package extensions when loaded;
# the device-agnostic KernelAbstractions kernels live in the GPU extension.
const _GPU_DEVICE = Ref{Any}(nothing)

"""
    _register_gpu_device!(device; supports_f64)

Called by the Metal/CUDA package extensions on load to register a
KernelAbstractions device for the `:gpu` backend.
"""
function _register_gpu_device!(device; supports_f64::Bool)
    _GPU_DEVICE[] = (device=device, f64=supports_f64)
    nothing
end

# Hooks specialized for the GPU backend in ff_gpu.jl (no-ops for CPU backends).
_backend_init!(::CompiledForceField) = nothing        # called at end of compile
_backend_upload_pairs!(::CompiledForceField) = nothing # called at end of rebuild_pairlist!

@inline _row(cff::CompiledForceField, atom) = cff.idx2row[atom.idx]

# --- bonded lowering (done once) --------------------------------------------
function _lower_bonded!(cff::CompiledForceField{T,Acc}) where {T,Acc}
    for c in cff.ff.components
        _lower_component!(cff, c)
    end
    cff
end

_lower_component!(::CompiledForceField, ::AbstractForceFieldComponent) = nothing

function _lower_component!(cff::CompiledForceField{T,Acc}, c::QuadraticStretchComponent{T}) where {T,Acc}
    s = cff.stretch
    for qbs in c.stretches
        push!(s.i, _row(cff, qbs.a1));
        push!(s.j, _row(cff, qbs.a2))
        push!(s.r0, Acc(qbs.r0));
        push!(s.k, Acc(qbs.k))
    end
    nothing
end

function _lower_component!(cff::CompiledForceField{T,Acc}, c::QuadraticBendComponent{T}) where {T,Acc}
    b = cff.bend
    for qab in c.bends
        push!(b.i, _row(cff, qab.a1));
        push!(b.j, _row(cff, qab.a2));
        push!(b.l, _row(cff, qab.a3))
        push!(b.θ0, Acc(qab.θ₀));
        push!(b.k, Acc(qab.k))
    end
    nothing
end

function _lower_component!(cff::CompiledForceField{T,Acc}, c::TorsionComponent{T}) where {T,Acc}
    _lower_torsions!(cff, cff.proper, c.proper_torsions)
    _lower_torsions!(cff, cff.improper, c.improper_torsions)
    nothing
end

function _lower_torsions!(cff::CompiledForceField{T,Acc}, t::TorsionArrays{Acc}, torsions) where {T,Acc}
    for ct in torsions
        push!(t.i, _row(cff, ct.a1));
        push!(t.j, _row(cff, ct.a2))
        push!(t.k, _row(cff, ct.a3));
        push!(t.l, _row(cff, ct.a4))
        for n in eachindex(ct.V)
            push!(t.V, Acc(ct.V[n]));
            push!(t.ϕ0, Acc(ct.ϕ₀[n]))
            push!(t.f, Int32(ct.f[n]));
            push!(t.div, Int32(ct.div[n]))
        end
        push!(t.term_offset, Int32(length(t.V) + 1))
    end
    nothing
end

# --- nonbonded lowering (rebuilt on demand) ---------------------------------
"""
    rebuild_pairlist!(cff::CompiledForceField)

Recompute the cached nonbonded pair list (vdW / hydrogen-bond / electrostatic)
from the force field's current coordinates. This is the expensive operation;
during minimization it is invoked only periodically (see `update_frequency` /
`max_displacement`).
"""
# Diagnostic: how many full pair-list rebuilds have happened. Rebuilds are serial
# and O(#pairs), so a high count during minimization explains poor thread-scaling
# (each `prepare_eval!` past the Verlet skin triggers one). Reset/read around a run.
const _REBUILD_COUNT = Ref(0)
reset_rebuild_count!() = (_REBUILD_COUNT[] = 0; nothing)
rebuild_count() = _REBUILD_COUNT[]

function rebuild_pairlist!(cff::CompiledForceField{T,Acc}) where {T,Acc}
    _REBUILD_COUNT[] += 1
    ff = cff.ff
    nbc = nothing
    for c in ff.components
        c isa NonBondedComponent && (nbc = c)
    end
    nbc === nothing && return cff

    # Push the evaluator's current coordinates into the canonical table so the
    # reference `update!` (which reads atoms(system).r and the CellListMap
    # neighbor search) sees the right positions. This also keeps the table
    # current at every rebuild (good for observables / visualization).
    sync_to_table!(cff; forces=false)

    # Re-run the reference nonbonded update (applies exclusions, HB/vdW
    # classification, 1-4 scaling, switching assignment). Single-sourced.
    update!(nbc)

    cff.vdw_sw = _convert_sw(Acc, nbc.cache.vdw_switching_function)
    cff.es_sw = _convert_sw(Acc, nbc.cache.es_switching_function)
    cff.distance_dependent_dielectric = ff.options[:distance_dependent_dielectric]::Bool

    _lower_lj!(cff, cff.lj, nbc.lj_interactions)
    _lower_lj!(cff, cff.hb, nbc.hydrogen_bonds)
    _lower_es!(cff, cff.es, nbc.electrostatic_interactions)

    # snapshot positions for the displacement check
    resize!(cff.r_at_last_build, cff.natoms)
    @inbounds for k in 1:cff.natoms
        cff.r_at_last_build[k] = cff.r[k]
    end
    cff.iters_since_build = 0
    _backend_upload_pairs!(cff)   # GPU extensions re-upload the new pair arrays
    cff
end

# largest atomic displacement since the last pair-list build
function _max_displacement(cff::CompiledForceField{T,Acc}) where {T,Acc}
    maxd2 = zero(Acc)
    @inbounds for k in 1:cff.natoms
        d2 = squared_norm(cff.r[k] - cff.r_at_last_build[k])
        d2 > maxd2 && (maxd2 = d2)
    end
    sqrt(maxd2)
end

"""
    maybe_rebuild!(cff::CompiledForceField) -> Bool

Rebuild the nonbonded pair list if atoms have drifted past the Verlet skin
(`max_displacement`) since the last build, or — when `update_frequency > 0` —
every `update_frequency` calls. Returns `true` if a rebuild happened. This is
the BALL-style decoupling of pair-list maintenance from per-evaluation cost.
"""
function maybe_rebuild!(cff::CompiledForceField{T,Acc}) where {T,Acc}
    need = _max_displacement(cff) > cff.max_displacement
    if cff.update_frequency > 0
        cff.iters_since_build += 1
        need |= cff.iters_since_build >= cff.update_frequency
    end
    need && rebuild_pairlist!(cff)
    need
end

"""
    prepare_eval!(cff, r)

Load the flat optimizer coordinate vector `r` into the evaluator and refresh the
pair list if necessary. Call before `compute_energy!`/`compute_forces!`.
"""
@inline function prepare_eval!(cff::CompiledForceField, r::AbstractVector)
    set_positions_flat!(cff, r)
    maybe_rebuild!(cff)
    cff
end

_convert_sw(::Type{Acc}, sw::CubicSwitchingFunction) where Acc =
    CubicSwitchingFunction{Acc}(Acc(sw.cutoff), Acc(sw.cuton))

function _lower_lj!(cff::CompiledForceField{T,Acc}, p::PairLJArrays{Acc}, interactions) where {T,Acc}
    empty!(p.i);
    empty!(p.j);
    empty!(p.A);
    empty!(p.B);
    empty!(p.scaling)
    for lji in interactions
        push!(p.i, _row(cff, lji.a1));
        push!(p.j, _row(cff, lji.a2))
        push!(p.A, Acc(lji.A));
        push!(p.B, Acc(lji.B));
        push!(p.scaling, Acc(lji.scaling_factor))
    end
    nothing
end

function _lower_es!(cff::CompiledForceField{T,Acc}, p::PairESArrays{Acc}, interactions) where {T,Acc}
    empty!(p.i);
    empty!(p.j);
    empty!(p.q1q2);
    empty!(p.scaling)
    for esi in interactions
        push!(p.i, _row(cff, esi.a1));
        push!(p.j, _row(cff, esi.a2))
        push!(p.q1q2, Acc(esi.q1q2));
        push!(p.scaling, Acc(esi.scaling_factor))
    end
    nothing
end

# ----------------------------------------------------------------------------
#  Coordinate / force synchronisation with the canonical atom table.
#  The hot loop never touches the table; these run only at sync points
#  (notify callbacks, end of optimization).
# ----------------------------------------------------------------------------
function sync_positions_from_table!(cff::CompiledForceField{T,Acc}) where {T,Acc}
    at = atoms(cff.ff.system)
    @inbounds for k in 1:cff.natoms
        cff.r[k] = Vector3{Acc}(at.r[k])
    end
    cff
end

"""
    sync_to_table!(cff; forces=true)

Write the evaluator's coordinates (and, by default, forces) back into the
canonical atom table so the rest of the library and any registered observables
see consistent state.
"""
function sync_to_table!(cff::CompiledForceField{T,Acc}; forces::Bool=true) where {T,Acc}
    at = atoms(cff.ff.system)
    @inbounds for k in 1:cff.natoms
        at.r[k] = Vector3{T}(cff.r[k])
    end
    if forces
        @inbounds for k in 1:cff.natoms
            at.F[k] = Vector3{T}(cff.F[k])
        end
    end
    cff
end

# pull positions from a flat optimizer vector (Float64) into the Acc SoA buffer
function set_positions_flat!(cff::CompiledForceField{T,Acc}, r::AbstractVector) where {T,Acc}
    @inbounds for k in 1:cff.natoms
        b = 3 * (k - 1)
        cff.r[k] = Vector3{Acc}(Acc(r[b+1]), Acc(r[b+2]), Acc(r[b+3]))
    end
    cff
end

# write the negative force (= energy gradient) into a flat optimizer vector
function gradient_flat!(grad::AbstractVector, cff::CompiledForceField{T,Acc}) where {T,Acc}
    @inbounds for k in 1:cff.natoms
        b = 3 * (k - 1)
        Fk = cff.F[k]
        grad[b+1] = -Fk[1];
        grad[b+2] = -Fk[2];
        grad[b+3] = -Fk[3]
    end
    grad
end


function _permute_rows!(v::AbstractVector, perm::AbstractVector{<:Integer})
    n = length(v)
    n == 0 && return v
    newv = similar(v)
    @inbounds for (dst, src) in enumerate(perm)
        newv[dst] = v[src]
    end
    copyto!(v, newv)
    v
end

function _permute_torsion_rows!(torsion::TorsionArrays{Acc}, perm::AbstractVector{<:Integer}) where {Acc}
    n = length(torsion.i)
    n == 0 && return torsion

    _permute_rows!(torsion.i, perm)
    _permute_rows!(torsion.j, perm)
    _permute_rows!(torsion.k, perm)
    _permute_rows!(torsion.l, perm)

    new_offsets = Vector{Int32}(undef, n + 1)
    new_offsets[1] = 1
    @inbounds for dst in 1:n
        src = Int(perm[dst])
        start = Int(torsion.term_offset[src])
        stop = Int(torsion.term_offset[src+1]) - 1
        new_offsets[dst+1] = new_offsets[dst] + Int32(stop - start + 1)
    end
    copyto!(torsion.term_offset, new_offsets)

    nterms = length(torsion.V)
    new_V = similar(torsion.V, nterms)
    new_ϕ0 = similar(torsion.ϕ0, nterms)
    new_f = similar(torsion.f, nterms)
    new_div = similar(torsion.div, nterms)

    pos = 1
    @inbounds for dst in 1:n
        src = Int(perm[dst])
        start = Int(torsion.term_offset[src])
        stop = Int(torsion.term_offset[src+1]) - 1
        len = stop - start + 1
        if len > 0
            copyto!(view(new_V, pos:(pos+len-1)), view(torsion.V, start:stop))
            copyto!(view(new_ϕ0, pos:(pos+len-1)), view(torsion.ϕ0, start:stop))
            copyto!(view(new_f, pos:(pos+len-1)), view(torsion.f, start:stop))
            copyto!(view(new_div, pos:(pos+len-1)), view(torsion.div, start:stop))
        end
        pos += len
    end

    copyto!(torsion.V, new_V)
    copyto!(torsion.ϕ0, new_ϕ0)
    copyto!(torsion.f, new_f)
    copyto!(torsion.div, new_div)

    torsion
end

function shuffle_interactions!(cff::CompiledForceField{T,Acc}; seed::Int=42) where {T,Acc}
    
    Random.seed!(seed)
    n = length(cff.stretch.i)
    if n > 0
        perm = randperm(n)
        _permute_rows!(cff.stretch.i, perm)
        _permute_rows!(cff.stretch.j, perm)
        _permute_rows!(cff.stretch.r0, perm)
        _permute_rows!(cff.stretch.k, perm)
    end

    n = length(cff.bend.i)
    if n > 0
        perm = randperm(n)
        _permute_rows!(cff.bend.i, perm)
        _permute_rows!(cff.bend.j, perm)
        _permute_rows!(cff.bend.l, perm)
        _permute_rows!(cff.bend.θ0, perm)
        _permute_rows!(cff.bend.k, perm)
    end

    n = length(cff.proper.i)
    if n > 0
        perm = randperm(n)
        _permute_torsion_rows!(cff.proper, perm)
    end

    n = length(cff.improper.i)
    if n > 0
        perm = randperm(n)
        _permute_torsion_rows!(cff.improper, perm)
    end

    n = length(cff.lj.i)
    if n > 0
        perm = randperm(n)
        _permute_rows!(cff.lj.i, perm)
        _permute_rows!(cff.lj.j, perm)
        _permute_rows!(cff.lj.A, perm)
        _permute_rows!(cff.lj.B, perm)
        _permute_rows!(cff.lj.scaling, perm)
    end

    n = length(cff.hb.i)
    if n > 0
        perm = randperm(n)
        _permute_rows!(cff.hb.i, perm)
        _permute_rows!(cff.hb.j, perm)
        _permute_rows!(cff.hb.A, perm)
        _permute_rows!(cff.hb.B, perm)
        _permute_rows!(cff.hb.scaling, perm)
    end

    n = length(cff.es.i)
    if n > 0
        perm = randperm(n)
        _permute_rows!(cff.es.i, perm)
        _permute_rows!(cff.es.j, perm)
        _permute_rows!(cff.es.q1q2, perm)
        _permute_rows!(cff.es.scaling, perm)
    end

    cff
end

function finalize_optimization!(cff::CompiledForceField{T,Acc}, u::AbstractVector) where {T,Acc}
    set_positions_flat!(cff, u)
    rebuild_pairlist!(cff)
    compute_forces!(cff)
    sync_to_table!(cff; forces=true)
    true
end
include("ff_kernels.jl")
