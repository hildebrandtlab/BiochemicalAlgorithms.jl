export Nucleotide, nucleotide_by_idx, nucleotides, nucleotides_df, eachnucleotide, nnucleotides,
    parent_nucleotide

"""
    $(TYPEDEF)

Mutable representation of an individual nucleotide in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
   Nucleotide(chain::Chain{T}, number::Int, name::String = "", properties::Properties = Properties())

Creates a new `Nucleotide{T}` in the given chain.
"""
@auto_hash_equals struct Nucleotide{T}
    _sys::System{T}
    _row::DataFrameRow
end

function Nucleotide(
    chain::Chain{T},
    number::Int,
    name::String = "",
    properties::Properties = Properties()
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(sys._nucleotides), (idx, number, name, properties, chain._row.molecule_id, chain.idx)
    nucleotide_by_idx(sys, idx)
end

function Base.getproperty(nuc::Nucleotide, name::Symbol)
    in(name, fieldnames(NucleotideTuple)) && return getproperty(getfield(nuc, :_row), name)
    getfield(nuc, name)
end

function Base.setproperty!(nuc::Nucleotide, name::Symbol, val)
    in(name, fieldnames(NucleotideTuple)) && return setproperty!(getfield(nuc, :_row), name, val)
    setfield!(nuc, name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", nuc::Nucleotide) = show(io, getfield(nuc, :_row))
@inline Base.show(io::IO, nuc::Nucleotide) = show(io, getfield(nuc, :_row))

@inline Base.parent(nuc::Nucleotide) = nuc._sys
@inline parent_system(nuc::Nucleotide) = parent(nuc)
@inline parent_molecule(nuc::Nucleotide) = molecule_by_idx(parent(nuc), nuc._row.molecule_id)
@inline parent_chain(nuc::Nucleotide) = chain_by_idx(parent(nuc), nuc._row.chain_id)

@doc raw"""
    parent_nucleotide(::Atom)

Returns the `Nucleotide{T}` containing the given atom. Returns `nothing` if no such nucleotide exists.
""" parent_nucleotide

"""
    $(TYPEDSIGNATURES)

Returns the `Nucleotide{T}` associated with the given `idx` in `sys`. Returns `nothing` if no such
nucleotide exists.
"""
@inline function nucleotide_by_idx(sys::System{T}, idx::Int) where T
    rn = _row_by_idx(sys._nucleotides, idx)
    isnothing(rn) ? nothing : Nucleotide{T}(sys, DataFrameRow(sys._nucleotides, rn, :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's nucleotides matching the given criteria. Fields
given as `nothing` are ignored. The returned `DataFrame` contains all public and private nucleotide fields.
"""
function _nucleotides(sys::System{T};
    molecule_id::Union{Nothing, Int} = nothing,
    chain_id::Union{Nothing, Int} = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return sys._nucleotides

    cols = Tuple{Symbol, Int}[]
    isnothing(molecule_id) || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)    || push!(cols, (:chain_id, chain_id))

    get(
        groupby(sys._nucleotides, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        DataFrame(_SystemNucleotideTuple[])
    )
end

"""
    nucleotides(::Chain)
    nucleotides(::Molecule)
    nucleotides(::Protein)
    nucleotides(::System)

Returns a `Vector{Nucleotide{T}}` containing all nucleotides of the given atom container.
"""
@inline function nucleotides(sys::System; kwargs...)
    collect(eachnucleotide(sys; kwargs...))
end

"""
    nucleotides_df(::Chain)
    nucleotides_df(::Molecule)
    nucleotides_df(::Protein)
    nucleotides_df(::System)

Returns a `SystemDataFrame{T}` containing all nucleotides of the given atom container.
"""
@inline function nucleotides_df(sys::System; kwargs...)
    SystemDataFrame(sys, view(_nucleotides(sys; kwargs...), :, 1:length(fieldnames(NucleotideTuple))))
end

"""
    eachnucleotide(::Chain)
    eachnucleotide(::Molecule)
    eachnucleotide(::Protein)
    eachnucleotide(::System)

Returns a `Nucleotide{T}` generator for all nucleotides of the given atom container.
"""
@inline function eachnucleotide(sys::System{T}; kwargs...) where T
    (Nucleotide{T}(sys, row) for row in eachrow(_nucleotides(sys; kwargs...)))
end

"""
    nnucleotides(::Chain)
    nnucleotides(::Molecule)
    nnucleotides(::Protein)
    nnucleotides(::System)

Returns the number of nucleotides in the given atom container.
"""
@inline function nnucleotides(sys::System; kwargs...)
    nrow(_nucleotides(sys; kwargs...))
end

#=
    Nucleotides
=#
@inline _nucleotides(mol::Molecule; kwargs...) = _nucleotides(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nucleotides(mol::Molecule; kwargs...) = nucleotides(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nucleotides_df(mol::Molecule; kwargs...) = nucleotides_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachnucleotide(mol::Molecule; kwargs...) = eachnucleotide(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nnucleotides(mol::Molecule; kwargs...) = nnucleotides(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain nucleotides
=#
@inline _nucleotides(chain::Chain; kwargs...) = _nucleotides(parent(chain); chain_id = chain.idx, kwargs...)
@inline nucleotides(chain::Chain; kwargs...) = nucleotides(parent(chain); chain_id = chain.idx, kwargs...)
@inline nucleotides_df(chain::Chain; kwargs...) = nucleotides_df(parent(chain); chain_id = chain.idx, kwargs...)
@inline eachnucleotide(chain::Chain; kwargs...) = eachnucleotide(parent(chain); chain_id = chain.idx, kwargs...)
@inline nnucleotides(chain::Chain; kwargs...) = nnucleotides(parent(chain); chain_id = chain.idx, kwargs...)

# FIXME currently not possible due to
# <https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/issues/26>
#
#@inline function Base.push!(chain::Chain{T}, nuc::NucleotideTuple) where T
#    push!(parent(chain)._nucleotides, (_with_idx(nuc, _next_idx(chain._sys))..., chain.idx))
#    chain
#end

#=
    Nucleotide atoms
=#
@inline _atoms(nuc::Nucleotide; kwargs...) = _atoms(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline atoms(nuc::Nucleotide; kwargs...) = atoms(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline atoms_df(nuc::Nucleotide; kwargs...) = atoms_df(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline eachatom(nuc::Nucleotide; kwargs...) = eachatom(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline natoms(nuc::Nucleotide; kwargs...) = natoms(parent(nuc); nucleotide_id = nuc.idx, kwargs...)

@inline function Base.push!(nuc::Nucleotide{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(parent(nuc), atom; molecule_id = nuc._row.molecule_id, chain_id = nuc._row.chain_id,
        nucleotide_id = nuc.idx, kwargs...)
    nuc
end

#=
    Nucleotide bonds
=#
@inline _bonds(nuc::Nucleotide; kwargs...) = _bonds(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline bonds(nuc::Nucleotide; kwargs...) = bonds(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline bonds_df(nuc::Nucleotide; kwargs...) = bonds_df(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline eachbond(nuc::Nucleotide; kwargs...) = eachbond(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline nbonds(nuc::Nucleotide; kwargs...) = nbonds(parent(nuc); nucleotide_id = nuc.idx, kwargs...)

@inline function Base.push!(nuc::Nucleotide, bond::Bond)
    push!(parent(nuc), bond)
    nuc
end

# TODO: we should come up with a better test than just checking the name
is_nucleotide(name::String) = name ∈ ["A", "C", "G", "T", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI"]
