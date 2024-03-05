export
    Fragment,
    FragmentTable,
    fragment_by_idx,
    fragments,
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
@inline Tables.rows(ft::FragmentTable) = ft

@inline function Tables.getcolumn(ft::FragmentTable, nm::Symbol)
    col = Tables.getcolumn(_fragments(ft), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _fragments(ft)._idx_map[idx], getfield(ft, :_idx))
    )
end

@inline Tables.getcolumn(ft::FragmentTable, i::Int) = Tables.getcolumn(ft, Tables.columnnames(ft)[i])
@inline Tables.columnnames(ft::FragmentTable) = Tables.columnnames(_fragments(ft))
@inline Tables.schema(ft::FragmentTable) = Tables.schema(_fragments(ft))

@inline Base.size(ft::FragmentTable) = (length(getfield(ft, :_idx)), length(_fragment_table_cols))
@inline Base.size(ft::FragmentTable, dim) = size(ft)[dim]
@inline Base.length(ft::FragmentTable) = size(ft, 1)

@inline function Base.getproperty(ft::FragmentTable, nm::Symbol)
    hasfield(typeof(ft), nm) && return getfield(ft, nm)
    Tables.getcolumn(ft, nm)
end

@inline function Base.setproperty!(ft::FragmentTable, nm::Symbol, val)
    if nm in _fragment_table_cols_priv || nm in _fragment_table_cols_set
        error("FragmentTable columns cannot be set directly! Did you mean to use broadcast assignment (.=)?")
    end
    if !hasfield(typeof(ft), nm)
        error("type FragmentTable has no field $nm")
    end
    setfield!(ft, nm, val)
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

@inline function Fragment(
    chain::Chain{T},
    number::Int;
    kwargs...
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(
        sys._fragments,
        _Fragment(number; idx = idx, kwargs...),
        chain.molecule_idx,
        chain.idx
    )
    fragment_by_idx(sys, idx)
end

@inline function Base.getproperty(frag::Fragment, name::Symbol)
    hasfield(typeof(frag), name) && return getfield(frag, name)
    getproperty(getfield(frag, :_row), name)
end

@inline function Base.setproperty!(frag::Fragment, name::Symbol, val)
    hasfield(typeof(frag), name) && return setfield!(frag, name, val)
    setproperty!(getfield(frag, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", frag::Fragment) = show(io, getfield(frag, :_row))
@inline Base.show(io::IO, frag::Fragment) = show(io, getfield(frag, :_row))

@inline Base.parent(frag::Fragment) = frag._sys
@inline parent_system(frag::Fragment) = parent(frag)
@inline parent_molecule(frag::Fragment) = molecule_by_idx(parent(frag), frag.molecule_idx)
@inline parent_chain(frag::Fragment) = chain_by_idx(parent(frag), frag.chain_idx)

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
    fragments(::System)

Returns a `FragmentTable{T}` containing all fragments of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to fragments belonging to the molecule with the given ID.
- `chain_idx::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to fragments belonging to the chain with the given ID.
"""
function fragments(sys::System{T};
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing
) where T
    isnothing(molecule_idx) && isnothing(chain_idx) && return FragmentTable{T}(sys, sys._fragments.idx)
    _filter_fragments(frag ->
        (isnothing(molecule_idx) || frag.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)    || frag.chain_idx    == something(chain_idx)),
        sys
    )
end

"""
    nfragments(::Chain)
    nfragments(::Molecule)
    nfragments(::System)

Returns the number of fragments in the given atom container.
"""
@inline function nfragments(sys::System; kwargs...)
    length(fragments(sys; kwargs...))
end

#=
    Molecule fragments
=#
@inline fragments(mol::Molecule; kwargs...) = fragments(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline nfragments(mol::Molecule; kwargs...) = nfragments(parent(mol); molecule_idx = mol.idx, kwargs...)

#=
    Chain fragments
=#
@inline fragments(chain::Chain; kwargs...) = fragments(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nfragments(chain::Chain; kwargs...) = nfragments(parent(chain); chain_idx = chain.idx, kwargs...)

"""
    push!(::Chain{T}, frag::Fragment{T})

Creates a copy of the given fragment in the given chain. The new fragment is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, frag::Fragment{T}) where T
    Fragment(chain, frag.number;
        name = frag.name,
        properties = frag.properties,
        flags = frag.flags
    )
    chain
end

#=
    Fragment atoms
=#
@inline atoms(frag::Fragment; kwargs...) = atoms(parent(frag); fragment_idx = frag.idx, kwargs...)
@inline natoms(frag::Fragment; kwargs...) = natoms(parent(frag); fragment_idx = frag.idx, kwargs...)

@inline function Atom(frag::Fragment, number::Int, element::ElementType; kwargs...)
    Atom(parent(frag), number, element;
        molecule_idx = frag.molecule_idx,
        chain_idx = frag.chain_idx,
        fragment_idx = frag.idx,
        kwargs...
    )
end

@inline function Base.push!(frag::Fragment{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(frag), atom;
        molecule_idx = frag.molecule_idx,
        chain_idx = frag.chain_idx,
        fragment_idx = frag.idx,
        kwargs...
    )
    frag
end

#=
    Fragment bonds
=#
@inline bonds(frag::Fragment; kwargs...) = bonds(parent(frag); fragment_idx = frag.idx, kwargs...)
@inline nbonds(frag::Fragment; kwargs...) = nbonds(parent(frag); fragment_idx = frag.idx, kwargs...)

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

            for f in fragments(c)
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

            for f in Iterators.reverse(fragments(c))
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

            for f in fragments(c)
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
