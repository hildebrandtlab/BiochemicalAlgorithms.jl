export
    Fragment,
    FragmentTable,
    fragment_by_idx,
    fragments,
    fragments_df,
    eachfragment,
    get_full_name,
    get_previous,
    get_next,
    is_3_prime,
    is_5_prime,
    is_amino_acid,
    is_c_terminal,
    is_n_terminal,
    is_nucleotide,
    nfragments,
    parent_fragment

@auto_hash_equals struct FragmentTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _fragments(ft::FragmentTable) = getproperty(getfield(ft, :_sys), :_fragments)

@inline Tables.istable(::Type{<: FragmentTable}) = true
@inline Tables.columnaccess(::Type{<: FragmentTable}) = true
@inline Tables.columns(ft::FragmentTable) = ft

@inline function Tables.getcolumn(ft::FragmentTable, nm::Symbol)
    col = Tables.getcolumn(_fragments(ft), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _fragments(ft)._idx_map[idx], getfield(ft, :_idx))
    )
end

@inline function Base.getproperty(ft::FragmentTable, nm::Symbol)
    hasfield(typeof(ft), nm) && return getfield(ft, nm)
    Tables.getcolumn(ft, nm)
end

@inline Tables.getcolumn(ft::FragmentTable, i::Int) = Tables.getcolumn(ft, Tables.columnnames(ft)[i])
@inline Tables.columnnames(ft::FragmentTable) = Tables.columnnames(_fragments(ft))
@inline Tables.schema(ft::FragmentTable) = Tables.schema(_fragments(ft))

@inline Base.size(ft::FragmentTable) = (length(getfield(ft, :_idx)), length(_fragment_table_cols))
@inline Base.size(ft::FragmentTable, dim) = size(ft)[dim]
@inline Base.length(ft::FragmentTable) = size(ft, 1)

function Base.push!(ft::FragmentTable, t::FragmentTuple, molecule_id::Int, chain_id::Int)
    sys = getfield(ft, :_sys)
    push!(sys._fragments, t, molecule_id, chain_id)
    push!(getfield(ft, :_idx), sys._curr_idx)
    ft
end

@inline function _filter_fragments(f::Function, sys::System{T}) where T
    FragmentTable{T}(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._fragments),
        :idx
    )))
end

@inline function Base.filter(f::Function, ft::FragmentTable)
    FragmentTable(getfield(ft, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, ft),
        :idx
    )))
end

@inline function Base.iterate(ft::FragmentTable, st = 1)
    st > length(ft) ?
        nothing :
        (fragment_by_idx(getfield(ft, :_sys), getfield(ft, :_idx)[st]), st + 1)
end
@inline Base.eltype(::FragmentTable{T}) where T = Fragment{T}
@inline Base.getindex(ft::FragmentTable{T}, i::Int) where T = fragment_by_idx(getfield(ft, :_sys), getfield(ft, :_idx)[i])
@inline Base.keys(ft::FragmentTable) = LinearIndices((length(ft),))

"""
    $(TYPEDEF)

Mutable representation of an individual fragment in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Fragment(
    chain::Chain{T},
    number::Int,
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Fragment{T}` in the given chain.
"""
@auto_hash_equals struct Fragment{T} <: AbstractAtomContainer{T}
    _sys::System{T}
    _row::_FragmentTableRow
end

function Fragment(
    chain::Chain{T},
    number::Int,
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(sys._fragments, FragmentTuple(number;
            idx = idx,
            name = name,
            properties = properties,
            flags = flags
        ), chain._row.molecule_id, chain.idx
    )
    fragment_by_idx(sys, idx)
end

@inline Tables.rows(ft::FragmentTable) = ft
@inline Tables.getcolumn(frag::Fragment, name::Symbol) = Tables.getcolumn(getfield(frag, :_row), name)

@inline function Base.getproperty(frag::Fragment, name::Symbol)
    hasfield(typeof(frag), name) && return getfield(frag, name)
    getproperty(getfield(frag, :_row), name)
end

@inline function Base.setproperty!(frag::Fragment, name::Symbol, val)
    hasfield(typeof(frag), name) && return setfield!(frag, name, val)
    setproperty!(getfield(frag, :_row), name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", frag::Fragment) = show(io, getfield(frag, :_row))
@inline Base.show(io::IO, frag::Fragment) = show(io, getfield(frag, :_row))

@inline Base.parent(frag::Fragment) = frag._sys
@inline parent_system(frag::Fragment) = parent(frag)
@inline parent_molecule(frag::Fragment) = molecule_by_idx(parent(frag), frag._row.molecule_id)
@inline parent_chain(frag::Fragment) = chain_by_idx(parent(frag), frag._row.chain_id)

@doc raw"""
    parent_fragment(::Atom)

Returns the `Fragment{T}` containing the given atom. Returns `nothing` if no such fragment exists.
""" parent_fragment

"""
    $(TYPEDSIGNATURES)

Returns the `Fragment{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
fragment exists.
"""
@inline function fragment_by_idx(sys::System{T}, idx::Int) where T
    Fragment{T}(sys, _row_by_idx(sys._fragments, idx))
end

"""
    fragments(::Chain)
    fragments(::Molecule)
    fragments(::Protein)
    fragments(::System)

Returns a `FragmentTable{T}` containing all fragments of the given atom container.

# Supported keyword arguments
 - `molecule_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to fragments belonging to the molecule with the given ID.
- `chain_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to fragments belonging to the chain with the given ID.
"""
function fragments(sys::System{T};
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return FragmentTable{T}(sys, sys._fragments.idx)
    _filter_fragments(frag ->
        (isnothing(molecule_id) || frag.molecule_id == something(molecule_id)) &&
        (isnothing(chain_id)    || frag.chain_id    == something(chain_id)),
        sys
    )
end

"""
    fragments_df(::Chain)
    fragments_df(::Molecule)
    fragments_df(::Protein)
    fragments_df(::System)

Returns a `DataFrame` containing all fragments of the given atom container.
"""
@inline function fragments_df(sys::System; kwargs...)
    DataFrame(fragments(sys; kwargs...))
end

"""
    eachfragment(::Chain)
    eachfragment(::Molecule)
    eachfragment(::Protein)
    eachfragment(::System)

Returns a `Fragment{T}` generator for all fragments of the given atom container.
"""
@inline function eachfragment(sys::System{T}; kwargs...) where T
    (frag for frag in fragments(sys; kwargs...))
end

"""
    nfragments(::Chain)
    nfragments(::Molecule)
    nfragments(::Protein)
    nfragments(::System)

Returns the number of fragments in the given atom container.
"""
@inline function nfragments(sys::System; kwargs...)
    length(fragments(sys; kwargs...))
end

#=
    Molecule fragments
=#
@inline fragments(mol::Molecule; kwargs...) = fragments(parent(mol); molecule_id = mol.idx, kwargs...)
@inline fragments_df(mol::Molecule; kwargs...) = fragments_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachfragment(mol::Molecule; kwargs...) = eachfragment(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nfragments(mol::Molecule; kwargs...) = nfragments(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain fragments
=#
@inline fragments(chain::Chain; kwargs...) = fragments(parent(chain); chain_id = chain.idx, kwargs...)
@inline fragments_df(chain::Chain; kwargs...) = fragments_df(parent(chain); chain_id = chain.idx, kwargs...)
@inline eachfragment(chain::Chain; kwargs...) = eachfragment(parent(chain); chain_id = chain.idx, kwargs...)
@inline nfragments(chain::Chain; kwargs...) = nfragments(parent(chain); chain_id = chain.idx, kwargs...)

"""
    push!(::Chain, frag::FragmentTuple)

Creates a new fragment in the given chain, based on the given tuple. The new fragment is automatically
assigned a new `idx`.
"""
@inline function Base.push!(chain::Chain, frag::FragmentTuple)
    sys = parent(chain)
    push!(sys._fragments,
        (; frag..., idx = _next_idx(sys)),
        chain._row.molecule_id,
        chain.idx
    )
    chain
end

#=
    Fragment atoms
=#
@inline atoms(frag::Fragment; kwargs...) = atoms(parent(frag); fragment_id = frag.idx, kwargs...)
@inline atoms_df(frag::Fragment; kwargs...) = atoms_df(parent(frag); fragment_id = frag.idx, kwargs...)
@inline eachatom(frag::Fragment; kwargs...) = eachatom(parent(frag); fragment_id = frag.idx, kwargs...)
@inline natoms(frag::Fragment; kwargs...) = natoms(parent(frag); fragment_id = frag.idx, kwargs...)

@inline function Base.push!(frag::Fragment{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(parent(frag), atom; molecule_id = frag._row.molecule_id, chain_id = frag._row.chain_id,
        fragment_id = frag.idx, kwargs...)
    frag
end

#=
    Fragment bonds
=#
@inline bonds(frag::Fragment; kwargs...) = bonds(parent(frag); fragment_id = frag.idx, kwargs...)
@inline bonds_df(frag::Fragment; kwargs...) = bonds_df(parent(frag); fragment_id = frag.idx, kwargs...)
@inline eachbond(frag::Fragment; kwargs...) = eachbond(parent(frag); fragment_id = frag.idx, kwargs...)
@inline nbonds(frag::Fragment; kwargs...) = nbonds(parent(frag); fragment_id = frag.idx, kwargs...)

@inline function is_amino_acid(frag::Fragment)
    is_amino_acid(frag.name)
end

@inline function is_nucleotide(frag::Fragment)
    is_nucleotide(frag.name)
end

# TODO: these should really be defined in Residue, not Fragment
@inline function is_n_terminal(frag::Fragment; precomputed=true)
    if (precomputed)
        has_flag(frag, :N_TERMINAL)
    else
        if is_amino_acid(frag)
            # find the first amino acid in the chain of this fragment
            c = parent_chain(frag)

            for f in eachfragment(c)
                if is_amino_acid(f)
                    return f == frag
                end
            end
        end

        false
    end
end

@inline function is_c_terminal(frag::Fragment; precomputed=true)
    if (precomputed)
        has_flag(frag, :C_TERMINAL)
    else
        if is_amino_acid(frag)
            # find the last amino acid in the chain of this fragment
            c = parent_chain(frag)

            for f in Iterators.reverse(eachfragment(c))
                if is_amino_acid(f)
                    return f == frag
                end
            end
        end

        false
    end
end

# TODO: these should really be defined in Nucleotide, not Fragment
@inline function is_3_prime(frag::Fragment; precomputed=true)
    if (precomputed)
        has_flag(frag, Symbol("3_PRIME"))
    else
        if is_nucleotide(frag)
            # find the first nucleotide in the chain of this fragment
            c = parent_chain(frag)

            for f in eachfragment(c)
                if is_nucleotide(f)
                    return f == frag
                end
            end
        end

        false
    end
end

@inline function is_5_prime(frag::Fragment; precomputed=true)
    if (precomputed)
        has_flag(frag, Symbol("5_PRIME"))
    else
        if is_nucleotide(frag)
            # find the last amino acid in the chain of this fragment
            c = parent_chain(frag)

            for f in reverse(fragments(c))
                if is_nucleotide(f)
                    return f == frag
                end
            end
        end

        false
    end
end

@inline function get_previous(frag::Fragment{T}) where {T<:Real}
    try
        prev_candidate = fragment_by_idx(parent(frag), frag.idx - 1)

        if parent_chain(prev_candidate) == parent_chain(frag)
            return prev_candidate
        end
    catch
    end

    nothing
end

@inline function get_next(frag::Fragment{T}) where {T<:Real}
    try
        prev_candidate = fragment_by_idx(parent(frag), frag.idx + 1)

        if parent_chain(prev_candidate) == parent_chain(frag)
            return prev_candidate
        end
    catch
    end

    nothing
end

function get_full_name(
        f::Fragment{T},
        type::FullNameType.T = FullNameType.ADD_VARIANT_EXTENSIONS) where {T<:Real}
    # retrieve the residue name and remove blanks
    full_name = strip(f.name)

    # if the variant extension should be added, do so
    if (   type == FullNameType.ADD_VARIANT_EXTENSIONS 
        || type == FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)
        suffix = "-"

        if is_n_terminal(f)
            suffix = "-N"
        end

        if is_c_terminal(f)
            suffix = "-C"
        end

        if (is_c_terminal(f) && is_n_terminal(f)) 
            suffix = "-M"
        end

        if (has_property(f, :PROPERTY__HAS_SSBOND))
            suffix *= "S"
        end

        if (suffix != "-")
            full_name *= suffix;
        end
    end

    if (   type == FullNameType.ADD_RESIDUE_ID 
        || type == FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)
  
        full_name *= string(f.number);
    end

    full_name
end
