export Atom, atoms, atoms_df, eachatom, natoms

struct Atom{T}
    sys::System{T}
    row::DataFrameRow
end

function Atom(
    sys::System{T},
    number::Int,
    element::ElementType,
    name::String = "",
    atomtype::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    has_velocity::Bool = false,
    has_force::Bool = false,
    properties::Properties = Properties();
    frame_id::Int = 1,
    molecule_id::MaybeInt = missing,
    chain_id::MaybeInt = missing,
    fragment_id::MaybeInt = missing,
    nucleotide_id::MaybeInt = missing,
    residue_id::MaybeInt = missing
) where T
    idx = _next_idx(sys)
    push!(sys.atoms, (idx, number, element, name, atomtype, r, v, F, has_velocity, has_force,
        properties, frame_id, molecule_id, chain_id, fragment_id, nucleotide_id, residue_id))
    _atom_by_idx(sys, idx)
end

function Atom(
    number::Int,
    element::ElementType,
    name::String = "",
    atomtype::String = "",
    r::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    v::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    F::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    has_velocity::Bool = false,
    has_force::Bool = false,
    properties::Properties = Properties();
    kwargs...
)
    Atom(default_system(), number, element, name, atomtype, r, v, F, has_velocity, has_force,
        properties; kwargs...)
end

@inline function Atom(sys::System{T}, atom::AtomTuple{T}; kwargs...) where T
    Atom(sys, ntuple(i -> atom[i+1], length(atom) - 1)...; kwargs...)
end

@inline function Atom(t::AtomTuple{Float32}; kwargs...)
    Atom(default_system(), t; kwargs...)
end

function Base.getproperty(atom::Atom, name::Symbol)
    in(name, fieldnames(AtomTuple)) && return getproperty(getfield(atom, :row), name)
    getfield(atom, name)
end

function Base.setproperty!(atom::Atom, name::Symbol, val)
    in(name, fieldnames(AtomTuple)) && return setproperty!(getfield(atom, :row), name, val)
    setfield!(atom, name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", atom::Atom) = show(io, getfield(atom, :row))
@inline Base.show(io::IO, atom::Atom) = show(io, getfield(atom, :row))

@inline Base.parent(atom::Atom) = atom.sys
@inline parent_system(atom::Atom) = parent(atom)

@inline function parent_molecule(atom::Atom) 
    ismissing(atom.row.molecule_id) ? nothing : _molecule_by_idx(atom.sys, atom.row.molecule_id)
end

@inline function parent_chain(atom::Atom)
    ismissing(atom.row.chain_id) ? nothing : _chain_by_idx(atom.sys, atom.row.chain_id)
end

@inline function parent_fragment(atom::Atom)
    ismissing(atom.row.fragment_id) ? nothing : _fragment_by_idx(atom.sys, atom.row.fragment_id)
end

@inline function parent_nucleotide(atom::Atom)
    ismissing(atom.row.nucleotide_id) ? nothing : _nucleotide_by_idx(atom.sys, atom.row.nucleotide_id)
end

@inline function parent_residue(atom::Atom)
    ismissing(atom.row.residue_id) ? nothing : _residue_by_idx(atom.sys, atom.row.residue_id)
end

@inline function _atom_by_idx(sys::System{T}, idx::Int) where T
    Atom{T}(sys, DataFrameRow(sys.atoms, findfirst(sys.atoms.idx .== idx), :))
end

function _atoms(sys::System{T};
    frame_id::Union{Nothing, Int} = 1,
    molecule_id::Union{Nothing, MaybeInt} = nothing,
    chain_id::Union{Nothing, MaybeInt} = nothing,
    fragment_id::Union{Nothing, MaybeInt} = nothing,
    nucleotide_id::Union{Nothing, MaybeInt} = nothing,
    residue_id::Union{Nothing, MaybeInt} = nothing
) where T
    cols = Tuple{Symbol, MaybeInt}[]
    isnothing(frame_id)      || push!(cols, (:frame_id, frame_id))
    isnothing(molecule_id)   || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)      || push!(cols, (:chain_id, chain_id))
    isnothing(fragment_id)   || push!(cols, (:fragment_id, fragment_id))
    isnothing(nucleotide_id) || push!(cols, (:nucleotide_id, nucleotide_id))
    isnothing(residue_id)    || push!(cols, (:residue_id, residue_id))

    get(
        groupby(sys.atoms, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        DataFrame(_SystemAtomTuple{T}[])
    )
end

@inline function atoms(sys::System; kwargs...)
    collect(eachatom(sys; kwargs...))
end

@inline function atoms_df(sys::System{T}; kwargs...) where T
    SystemDataFrame{T}(sys, view(_atoms(sys; kwargs...), :, 1:length(fieldnames(AtomTuple{T}))))
end

@inline function eachatom(sys::System{T}; kwargs...) where T
    (Atom{T}(sys, row) for row in eachrow(_atoms(sys; kwargs...)))
end

@inline function natoms(sys::System; kwargs...)
    nrow(_atoms(sys; kwargs...))
end

function Base.push!(sys::System{T}, atom::AtomTuple{T};
    frame_id::Int = 1,
    molecule_id::MaybeInt = missing,
    chain_id::MaybeInt = missing,
    fragment_id::MaybeInt = missing,
    nucleotide_id::MaybeInt = missing,
    residue_id::MaybeInt = missing
) where T
    push!(sys.atoms, 
        (_with_idx(atom, _next_idx(sys))...,
            frame_id, molecule_id, chain_id, fragment_id, nucleotide_id, residue_id)
    )
    sys
end
