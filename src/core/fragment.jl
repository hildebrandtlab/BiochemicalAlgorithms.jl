export Fragment, fragment_by_idx, fragments, fragments_df, eachfragment, nfragments, parent_fragment, 
    is_c_terminal, is_n_terminal, is_amino_acid, is_nucleotide, is_3_prime, is_5_prime,
    is_previous, get_previous, is_next, get_next, get_full_name

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
    _row::DataFrameRow
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
    push!(sys._fragments, (idx, number, name, properties, flags, chain._row.molecule_id, chain.idx))
    fragment_by_idx(sys, idx)
end

function Base.getproperty(frag::Fragment, name::Symbol)
    gp = () -> getproperty(getfield(frag, :_row), name)
    name === :idx        && return gp()::Int
    name === :number     && return gp()::Int
    name === :name       && return gp()::String
    name === :properties && return gp()::Properties
    name === :flags      && return gp()::Flags
    getfield(frag, name)
end

function Base.setproperty!(frag::Fragment, name::Symbol, val)
    in(name, fieldnames(FragmentTuple)) && return setproperty!(getfield(frag, :_row), name, val)
    setfield!(frag, name, val)
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
    Fragment{T}(sys, DataFrameRow(sys._fragments.df, _row_by_idx(sys._fragments, idx), :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's fragments matching the given criteria. Fields
given as `nothing` are ignored. The returned `DataFrame` contains all public and private fragment fields.
"""
function _fragments(sys::System{T};
        molecule_id::MaybeInt = nothing,
        chain_id::MaybeInt = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return sys._fragments.df

    cols = Tuple{Symbol, Int}[]
    isnothing(molecule_id) || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)    || push!(cols, (:chain_id, chain_id))

    get(
        groupby(sys._fragments.df, getindex.(cols, 1)),
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

Returns a `SubDataFrame` containing all fragments of the given atom container.
"""
@inline function fragments_df(sys::System; kwargs...)
    view(_fragments(sys; kwargs...), :, 1:length(fieldnames(FragmentTuple)))
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
@inline _fragments(mol::Molecule; kwargs...) = _fragments(parent(mol); molecule_id = mol.idx, kwargs...)
@inline fragments(mol::Molecule; kwargs...) = fragments(parent(mol); molecule_id = mol.idx, kwargs...)
@inline fragments_df(mol::Molecule; kwargs...) = fragments_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachfragment(mol::Molecule; kwargs...) = eachfragment(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nfragments(mol::Molecule; kwargs...) = nfragments(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain fragments
=#
@inline _fragments(chain::Chain; kwargs...) = _fragments(parent(chain); chain_id = chain.idx, kwargs...)
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
    push!(sys._fragments, (; frag..., 
        idx = _next_idx(sys),
        molecule_id = chain._row.molecule_id,
        chain_id = chain.idx
    ))
    chain
end

#=
    Fragment atoms
=#
@inline _atoms(frag::Fragment; kwargs...) = _atoms(parent(frag); fragment_id = frag.idx, kwargs...)
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
@inline _bonds(frag::Fragment; kwargs...) = _bonds(parent(frag); fragment_id = frag.idx, kwargs...)
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

@inline function is_previous(f1::Fragment{T}, f2::Fragment{T}) where {T<:Real}
    (parent_chain(f1) == parent_chain(f2)) && (f1.idx == f2.idx - 1)
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


@inline function is_next(f1::Fragment{T}, f2::Fragment{T}) where {T<:Real}
    (parent_chain(f1) == parent_chain(f2)) && (f1.idx == f2.idx + 1)
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
