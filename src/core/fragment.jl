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

"""
    Fragment{T} <: AbstractAtomContainer{T}

Mutable representation of an individual fragment in a system.

# Public fields
 - `idx::Int`
 - `number::Int`
 - `name::String`

# Private fields
 - `variant::FragmentVariantType`
 - `properties::Properties`
 - `flags::Flags`
 - `molecule_idx::Int`
 - `chain_idx::Int`

# Constructors
```julia
Fragment(
    chain::Chain{T},
    number::Int;
    # keyword arguments
    name::String = "",
    variant::FragmentVariantType = FragmentVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Fragment{T}` in the given chain.
"""
const Fragment{T} = AtomContainer{T, _FragmentTableRow}

@inline function Fragment(
    chain::Chain{T},
    number::Int;
    kwargs...
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(sys._fragments, idx, number, chain.molecule_idx, chain.idx; kwargs...)
    fragment_by_idx(sys, idx)
end

"""
    FragmentTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system fragments (or a subset thereof). Fragment tables can be
generated using [`fragments`](@ref) or filtered from other fragment tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `number::AbstractVector{Int}`
 - `name::AbstractVector{String}`

# Private columns
 - `variant::AbstractVector{FragmentVariantType}`
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `molecule_idx::AbstractVector{Int}`
 - `chain_idx::AbstractVector{Int}`
"""
const FragmentTable{T} = SystemComponentTable{T, Fragment{T}}

@inline function _filter_fragments(f::Function, sys::System{T}) where T
    FragmentTable{T}(sys, _filter_idx(f, sys._fragments))
end

@inline _table(sys::System{T}, ::Type{Fragment{T}}) where T = sys._fragments

@inline function _hascolumn(::Type{<: Fragment}, nm::Symbol)
    nm in _fragment_table_cols_set || nm in _fragment_table_cols_priv
end

@inline parent_molecule(frag::Fragment) = molecule_by_idx(parent(frag), frag.molecule_idx)
@inline parent_chain(frag::Fragment) = chain_by_idx(parent(frag), frag.chain_idx)

@doc raw"""
    parent_fragment(::Atom)

Returns the `Fragment{T}` containing the given atom. Returns `nothing` if no such fragment exists.
""" parent_fragment

"""
    fragment_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Fragment{T}

Returns the `Fragment{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
fragment exists.
"""
@inline function fragment_by_idx(sys::System{T}, idx::Int) where T
    Fragment{T}(sys, _row_by_idx(sys._fragments, idx))
end

@inline function fragment_by_idx(idx::Int)
    fragment_by_idx(default_system(), idx)
end

"""
    fragments(::Chain)
    fragments(::Molecule)
    fragments(::System = default_system())

Returns a `FragmentTable{T}` containing all fragments of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`
 - `chain_idx::MaybeInt = nothing`
All keyword arguments limit the results to fragments matching the given IDs. Keyword arguments set to
`nothing` are ignored.
"""
function fragments(sys::System{T} = default_system();
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
    nfragments(::System = default_system())

Returns the number of fragments in the given atom container.

# Supported keyword arguments
See [`fragments`](@ref)
"""
@inline function nfragments(sys::System = default_system(); kwargs...)
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
    push!(::Chain{T}, ::Fragment{T})

Creates a copy of the given fragment in the given chain. The new fragment is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, frag::Fragment{T}) where T
    Fragment(chain, frag.number;
        name = frag.name,
        variant = frag.variant,
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
