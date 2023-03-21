using AutoHashEquals
export Fragment, fragments, fragments_df, eachfragment, nfragments, parent_fragment, is_c_terminal, is_n_terminal, is_amino_acid

"""
    $(TYPEDEF)

Mutable representation of an individual fragment in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
    Fragment(chain::Chain{T}, number::Int, name::String = "", properties::Properties = Properties())

Creates a new `Fragment{T}` in the given chain.
"""
@auto_hash_equals struct Fragment{T}
    sys::System{T}
    row::DataFrameRow
end

function Fragment(
    chain::Chain{T}, 
    number::Int, 
    name::String = "", 
    properties::Properties = Properties()
) where T
    sys = chain.sys
    idx = _next_idx(sys)
    push!(sys.fragments, (idx, number, name, properties, chain.row.molecule_id, chain.idx))
    _fragment_by_idx(sys, idx)
end

function Base.getproperty(frag::Fragment, name::Symbol)
    in(name, fieldnames(FragmentTuple)) && return getproperty(getfield(frag, :row), name)
    getfield(frag, name)
end

function Base.setproperty!(frag::Fragment, name::Symbol, val)
    in(name, fieldnames(FragmentTuple)) && return setproperty!(getfield(frag, :row), name, val)
    setfield!(frag, name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", frag::Fragment) = show(io, getfield(frag, :row))
@inline Base.show(io::IO, frag::Fragment) = show(io, getfield(frag, :row))

@inline Base.parent(frag::Fragment) = frag.sys
@inline parent_system(frag::Fragment) = parent(frag)
@inline parent_molecule(frag::Fragment) = _molecule_by_idx(frag.sys, frag.row.molecule_id)
@inline parent_chain(frag::Fragment) = _chain_by_idx(frag.sys, frag.row.chain_id)

@doc raw"""
    parent_fragment(::Atom)

Returns the `Fragment{T}` containing the given atom. Returns `nothing` if no such fragment exists.
""" parent_fragment

"""
    $(TYPEDSIGNATURES)

Returns the `Fragment{T}` associated with the given `idx` in `sys`.
"""
@inline function _fragment_by_idx(sys::System{T}, idx::Int) where T
    Fragment{T}(sys, DataFrameRow(sys.fragments, findfirst(sys.fragments.idx .== idx), :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's fragments matching the given criteria. Fields
given as `nothing` are ignored. The returned `DataFrame` contains all public and private fragment fields.
"""
function _fragments(sys::System{T};
        molecule_id::Union{Nothing, Int} = nothing,
        chain_id::Union{Nothing, Int} = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return sys.fragments

    cols = Tuple{Symbol, Int}[]
    isnothing(molecule_id) || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)    || push!(cols, (:chain_id, chain_id))

    get(
        groupby(sys.fragments, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        DataFrame(_SystemFragmentTuple[])
    )
end

"""
    fragments(::Chain)
    fragments(::Molecule)
    fragments(::Protein)
    fragments(::System)

Returns a `Vector{Fragment{T}}` containing all fragments of the given atom container.
"""
@inline function fragments(sys::System; kwargs...)
    collect(eachfragment(sys; kwargs...))
end

"""
    fragments_df(::Chain)
    fragments_df(::Molecule)
    fragments_df(::Protein)
    fragments_df(::System)

Returns a `SystemDataFrame{T}` containing all fragments of the given atom container.
"""
@inline function fragments_df(sys::System; kwargs...)
    SystemDataFrame(sys, view(_fragments(sys; kwargs...), :, 1:length(fieldnames(FragmentTuple))))
end

"""
    eachfragment(::Chain)
    eachfragment(::Molecule)
    eachfragment(::Protein)
    eachfragment(::System)

Returns a `Fragment{T}` generator for all fragments of the given atom container.
"""
@inline function eachfragment(sys::System{T}; kwargs...) where T
    (Fragment{T}(sys, row) for row in eachrow(_fragments(sys; kwargs...)))
end

"""
    nfragments(::Chain)
    nfragments(::Molecule)
    nfragments(::Protein)
    nfragments(::System)

Returns the number of fragments in the given atom container.
"""
@inline function nfragments(sys::System; kwargs...)
    nrow(_fragments(sys; kwargs...))
end

#=
    Molecule fragments
=#
@inline _fragments(mol::Molecule; kwargs...) = _fragments(mol.sys; molecule_id = mol.idx, kwargs...)
@inline fragments(mol::Molecule; kwargs...) = fragments(mol.sys; molecule_id = mol.idx, kwargs...)
@inline fragments_df(mol::Molecule; kwargs...) = fragments_df(mol.sys; molecule_id = mol.idx, kwargs...)
@inline eachfragment(mol::Molecule; kwargs...) = eachfragment(mol.sys; molecule_id = mol.idx, kwargs...)
@inline nfragments(mol::Molecule; kwargs...) = nfragments(mol.sys; molecule_id = mol.idx, kwargs...)

#=
    Chain fragments
=#
@inline _fragments(chain::Chain; kwargs...) = _fragments(chain.sys; chain_id = chain.idx, kwargs...)
@inline fragments(chain::Chain; kwargs...) = fragments(chain.sys; chain_id = chain.idx, kwargs...)
@inline fragments_df(chain::Chain; kwargs...) = fragments_df(chain.sys; chain_id = chain.idx, kwargs...)
@inline eachfragment(chain::Chain; kwargs...) = eachfragment(chain.sys; chain_id = chain.idx, kwargs...)
@inline nfragments(chain::Chain; kwargs...) = nfragments(chain.sys; chain_id = chain.idx, kwargs...)

"""
    push!(::Chain, frag::FragmentTuple)

Creates a new fragment in the given chain, based on the given tuple. The new fragment is automatically
assigned a new `idx`.
"""
@inline function Base.push!(chain::Chain, frag::FragmentTuple)
    push!(chain.sys.fragments, (_with_idx(frag, _next_idx(chain.sys))..., chain.idx))
    chain
end

#=
    Fragment atoms
=#
@inline _atoms(frag::Fragment; kwargs...) = _atoms(frag.sys; fragment_id = frag.idx, kwargs...)
@inline atoms(frag::Fragment; kwargs...) = atoms(frag.sys; fragment_id = frag.idx, kwargs...)
@inline atoms_df(frag::Fragment; kwargs...) = atoms_df(frag.sys; fragment_id = frag.idx, kwargs...)
@inline eachatom(frag::Fragment; kwargs...) = eachatom(frag.sys; fragment_id = frag.idx, kwargs...)
@inline natoms(frag::Fragment; kwargs...) = natoms(frag.sys; fragment_id = frag.idx, kwargs...)

@inline function Base.push!(frag::Fragment{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(frag.sys, atom; molecule_id = frag.row.molecule_id, chain_id = frag.row.chain_id,
        fragment_id = frag.idx, kwargs...)
    frag
end

#=
    Fragment bonds
=#
@inline _bonds(frag::Fragment; kwargs...) = _bonds(frag.sys; fragment_id = frag.idx, kwargs...)
@inline bonds(frag::Fragment; kwargs...) = bonds(frag.sys; fragment_id = frag.idx, kwargs...)
@inline bonds_df(frag::Fragment; kwargs...) = bonds_df(frag.sys; fragment_id = frag.idx, kwargs...)
@inline eachbond(frag::Fragment; kwargs...) = eachbond(frag.sys; fragment_id = frag.idx, kwargs...)
@inline nbonds(frag::Fragment; kwargs...) = nbonds(frag.sys; fragment_id = frag.idx, kwargs...)

@inline function Base.push!(frag::Fragment, bond::Bond)
    push!(frag.sys, bond)
    frag
end

@inline function is_amino_acid(frag::Fragment)
    is_amino_acid(frag.name)
end

# TODO: these should really be defined in Residue, not Fragment

@inline function is_n_terminal(frag::Fragment)
    if is_amino_acid(frag)
        # find the first amino acid in the chain of this fragment
        c = parent_chain(frag)

        for f in fragments(c)
            if is_amino_acid(f)
                return f == frag
            end
        end
    end

    false
end

@inline function is_c_terminal(frag::Fragment)
    if is_amino_acid(frag)
        # find the last amino acid in the chain of this fragment
        c = parent_chain(frag)

        for f in reverse(fragments(c))
            if is_amino_acid(f)
                return f == frag
            end
        end
    end

    false
end