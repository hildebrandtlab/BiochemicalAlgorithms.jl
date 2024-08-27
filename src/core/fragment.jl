export
    Fragment,
    FragmentTable,
    Nucleotide,
    Residue,
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
    isnucleotide,
    isresidue,
    nfragments,
    nnucleotides,
    nresidues,
    nucleotide_by_idx,
    nucleotides,
    parent_fragment,
    parent_nucleotide,
    parent_residue,
    residue_by_idx,
    residues

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
 - `secondary_structure_idx::Int`

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

```julia
Fragment(
    secondary_structure::SecondaryStructure{T},
    number::Int;
    # keyword arguments
    name::String = "",
    variant::FragmentVariantType = FragmentVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Fragment{T}` in the given secondary structure.

"""
const Fragment{T} = AtomContainer{T, :Fragment}

@inline function Fragment(
    chain::Chain{T},
    number::Int;
    kwargs...
) where T
    sys = parent(chain)
    idx = _next_idx!(sys)
    push!(sys._fragments, idx, number, chain.molecule_idx, chain.idx; kwargs...)
    fragment_by_idx(sys, idx)
end

@inline function Fragment(
    ss::SecondaryStructure{T},
    number::Int;
    kwargs...
) where T
    Fragment(parent_chain(ss), number; secondary_structure_idx = ss.idx, kwargs...)
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
    _hascolumn(_FragmentTable, nm)
end

@inline parent_molecule(frag::Fragment) = molecule_by_idx(parent(frag), frag.molecule_idx)
@inline parent_chain(frag::Fragment) = chain_by_idx(parent(frag), frag.chain_idx)
@inline parent_secondary_structure(frag::Fragment) = secondary_structure_by_idx(parent(frag), frag.secondary_structure_idx)

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
    _rowno_by_idx(_table(sys, Fragment{T}), idx) # check idx
    Fragment{T}(sys, idx)
end

@inline function fragment_by_idx(idx::Int)
    fragment_by_idx(default_system(), idx)
end

"""
    fragments(::Chain)
    fragments(::SecondaryStructure)
    fragments(::Molecule)
    fragments(::System = default_system())

Returns a `FragmentTable{T}` containing all fragments of the given atom container.

# Supported keyword arguments
 - `variant::Union{Nothing, FragmentVariantType} = nothing`
 - `molecule_idx::MaybeInt = nothing`
 - `chain_idx::MaybeInt = nothing`
 - `secondary_structure_idx::MaybeInt = nothing`
All keyword arguments limit the results to fragments matching the given IDs or variant type.
Keyword arguments set to `nothing` are ignored.
"""
function fragments(sys::System{T} = default_system();
    variant::Union{Nothing, FragmentVariantType} = nothing,
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing,
    secondary_structure_idx::MaybeInt = nothing
) where T
    isnothing(variant) &&
        isnothing(molecule_idx) &&
        isnothing(chain_idx) &&
        isnothing(secondary_structure_idx) &&
        return FragmentTable{T}(sys, copy(sys._fragments.idx))
    _filter_fragments(frag ->
        (isnothing(variant)                 || frag.variant      == something(variant)) &&
        (isnothing(molecule_idx)            || frag.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)               || frag.chain_idx    == something(chain_idx)) &&
        (isnothing(secondary_structure_idx) || frag.secondary_structure_idx == something(secondary_structure_idx)),
        sys
    )
end

"""
    nfragments(::Chain)
    nfragments(::SecondaryStructure)
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

"""
    fragments(::SecondaryStructureTable)
    fragments(::ChainTable)
    fragments(::MoleculeTable)

Returns a `FragmentTable{T}` containing all fragments of the given table.

# Supported keyword arguments
 - `variant::Union{Nothing, FragmentVariantType} = nothing`
   Any value other than `nothing` limits the results to the matching variant type.
"""
@inline function fragments(mt::MoleculeTable;
    variant::Union{Nothing, FragmentVariantType} = nothing
)
    idx = Set(mt.idx)
    isnothing(variant) ?
        _filter_fragments(frag -> frag.molecule_idx in idx, mt._sys) :
        _filter_fragments(frag -> frag.molecule_idx in idx && frag.variant == something(variant), mt._sys)
end

"""
    nfragments(::ChainTable)
    nfragments(::FragmentTable)
    nfragments(::MoleculeTable)

Returns the number of fragments belonging to the given table.

# Supported keyword arguments
See [`fragments`](@ref fragments)
"""
@inline function nfragments(
    ft::FragmentTable;
    variant::Union{Nothing, FragmentVariantType} = nothing
)
    isnothing(variant) ?
        length(ft) :
        length(filter(frag -> frag.variant == something(variant), ft))
end

@inline function nfragments(mt::MoleculeTable; kwargs...)
    length(fragments(mt; kwargs...))
end

#=
    Chain fragments
=#
@inline fragments(chain::Chain; kwargs...) = fragments(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nfragments(chain::Chain; kwargs...) = nfragments(parent(chain); chain_idx = chain.idx, kwargs...)

@inline function fragments(ct::ChainTable;
    variant::Union{Nothing, FragmentVariantType} = nothing
)
    idx = Set(ct.idx)
    isnothing(variant) ?
        _filter_fragments(frag -> frag.chain_idx in idx, ct._sys) :
        _filter_fragments(frag -> frag.chain_idx in idx && frag.variant == something(variant), ct._sys)
end

@inline nfragments(ct::ChainTable; kwargs...) = length(fragments(ct; kwargs...))

#=
    SecondaryStructure fragments
=#
@inline fragments(ss::SecondaryStructure; kwargs...) = fragments(parent(ss); secondary_structure_idx = ss.idx, kwargs...)
@inline nfragments(ss::SecondaryStructure; kwargs...) = nfragments(parent(ss); secondary_structure_idx = ss.idx, kwargs...)

@inline function fragments(st::SecondaryStructureTable;
        variant::Union{Nothing, FragmentVariantType} = nothing
)
    idx = Set(st.idx)
    isnothing(variant) ?
        _filter_fragments(frag -> frag.secondary_structure_idx in idx, st._sys) :
        _filter_fragments(frag -> frag.secondary_structure_idx in idx && frag.variant == something(variant), st._sys)
end

@inline nfragments(st::SecondaryStructureTable; kwargs...) = length(fragments(st; kwargs...))

"""
    delete!(::Fragment)
    delete!(::FragmentTable)
    delete!(::FragmentTable, idx::Int)

Removes the given fragment(s) from the associated system.

# Supported keyword arguments
 - `keep_atoms::Bool = false`
   Determines whether associated atoms (and their bonds) are removed as well
"""
function Base.delete!(frag::Fragment; keep_atoms::Bool = false)
    keep_atoms ? atoms(frag).fragment_idx .= Ref(nothing) : delete!(atoms(frag))
    delete!(parent(frag)._fragments, frag.idx)
    nothing
end

function Base.delete!(ft::FragmentTable; keep_atoms::Bool = false)
    keep_atoms ? atoms(ft).fragment_idx .= Ref(nothing) : delete!(atoms(ft))
    delete!(_table(ft), ft._idx)
    empty!(ft._idx)
    ft
end

function Base.delete!(ft::FragmentTable, idx::Int; kwargs...)
    idx in ft._idx || throw(KeyError(idx))
    delete!(fragment_by_idx(ft._sys, idx); kwargs...)
    deleteat!(ft._idx, findall(i -> i == idx, ft._idx))
    ft
end

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

"""
    push!(::SecondaryStructure{T}, ::Fragment{T})

Creates a copy of the given fragment in the given SecondaryStructure. The new fragment is automatically assigned a
new `idx`.
"""
@inline function Base.push!(ss::SecondaryStructure{T}, frag::Fragment{T}) where T
    Fragment(ss, frag.number;
        name = frag.name,
        variant = frag.variant,
        properties = frag.properties,
        flags = frag.flags
    )
    ss
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

@inline function atoms(ft::FragmentTable)
    idx = Set(ft._idx)
    _filter_atoms(atom -> atom.fragment_idx in idx, ft._sys)
end

@inline natoms(ft::FragmentTable) = length(atoms(ft))

#=
    Fragment bonds
=#
@inline bonds(frag::Fragment; kwargs...) = bonds(parent(frag); fragment_idx = frag.idx, kwargs...)
@inline nbonds(frag::Fragment; kwargs...) = nbonds(parent(frag); fragment_idx = frag.idx, kwargs...)

@inline bonds(ft::FragmentTable) = bonds(atoms(ft))
@inline nbonds(ft::FragmentTable) = nbonds(atoms(ft))

#=
    Variant: Nucleotide
=#
"""
    Nucleotide(chain::Chain, number::Int)

`Fragment` constructor defaulting to the [`FragmentVariant.Nucleotide`](@ref FragmentVariant) variant.

# Supported keyword arguments
See [`Fragment`](@ref)
"""
@inline function Nucleotide(chain::Chain, number::Int; kwargs...)
    Fragment(chain, number; variant = FragmentVariant.Nucleotide, kwargs...)
end


"""
    Nucleotide(ss::SecondaryStructure, number::Int)

`Fragment` constructor defaulting to the [`FragmentVariant.Nucleotide`](@ref FragmentVariant) variant.

# Supported keyword arguments
See [`Fragment`](@ref)
"""
@inline function Nucleotide(ss::SecondaryStructure, number::Int; kwargs...)
    Fragment(ss, number; variant = FragmentVariant.Nucleotide, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Returns `true` if the given fragment is a [`FragmentVariant.Nucleotide`](@ref FragmentVariant),
`false` otherwise.
"""
@inline function isnucleotide(frag::Fragment)
    frag.variant === FragmentVariant.Nucleotide
end

"""
    nucleotide_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Fragment{T}

Returns the `Fragment{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
fragment exists or if the fragment is not a [`FragmentVariant.Nucleotide`](@ref FragmentVariant).
"""
@inline function nucleotide_by_idx(sys::System, idx::Int)
    frag = fragment_by_idx(sys, idx)
    isnucleotide(frag) || throw(KeyError(idx))
    frag
end

@inline function nucleotide_by_idx(idx::Int)
    nucleotide_by_idx(default_system(), idx)
end

"""
    nucleotides(::Chain)
    nucleotides(::ChainTable)
    nucleotides(::SecondaryStructure)
    nucleotides(::SecondaryStructureTable)
    nucleotides(::Molecule)
    nucleotides(::MoleculeTable)
    nucleotides(::System = default_system())

Returns a `FragmentTable{T}` containing all [`FragmentVariant.Nucleotide`](@ref FragmentVariant)
fragments of the given atom container or table.

# Supported keyword arguments
See [`fragments`](@ref)
"""
@inline function nucleotides(
    ac::Union{Chain, ChainTable, SecondaryStructure, SecondaryStructureTable, Molecule, MoleculeTable, System} = default_system();
    kwargs...
)
    fragments(ac; variant = FragmentVariant.Nucleotide, kwargs...)
end

"""
    nnucleotides(::Chain)
    nnucleotides(::SecondaryStructure)
    nnucleotides(::Molecule)
    nnucleotides(::System = default_system())

Returns the number of [`FragmentVariant.Nucleotide`](@ref FragmentVariant) fragments in the given
atom container.

# Supported keyword arguments
See [`fragments`](@ref)
"""
@inline function nnucleotides(
    ac::Union{Chain, ChainTable, SecondaryStructure, SecondaryStructureTable, Molecule, MoleculeTable, System} = default_system();
    kwargs...
)
    nfragments(ac; variant = FragmentVariant.Nucleotide, kwargs...)
end

"""
    nnucleotides(::ChainTable)
    nnucleotides(::SecondaryStructureTable)
    nnucleotides(::FragmentTable)
    nnucleotides(::MoleculeTable)

Returns the number of [`FragmentVariant.Nucleotide`](@ref FragmentVariant) fragments in the given
table.
"""
@inline function nnucleotides(ft::FragmentTable)
    length(filter(isnucleotide, ft))
end

@inline function nnucleotides(ct::Union{ChainTable, SecondaryStructureTable, MoleculeTable})
    nfragments(ct; variant = FragmentVariant.Nucleotide)
end

"""
    parent_nucleotide(::Atom)

Returns the `Fragment{T}` containing the given atom. Returns `nothing` if no such fragment exists or
if the fragment is not a [`FragmentVariant.Nucleotide`](@ref FragmentVariant).
"""
@inline function parent_nucleotide(atom::Atom)
    frag = parent_fragment(atom)
    isnothing(frag) || !isnucleotide(frag) ? nothing : frag
end

#=
    Variant: Residue
=#
"""
    Residue(chain::Chain, number::Int)

`Fragment` constructor defaulting to the [`FragmentVariant.Residue`](@ref FragmentVariant) variant.

# Supported keyword arguments
See [`Fragment`](@ref)
"""
@inline function Residue(chain::Chain, number::Int; kwargs...)
    Fragment(chain, number; variant = FragmentVariant.Residue, kwargs...)
end

"""
    Residue(ss::SecondaryStructure, number::Int)

`Fragment` constructor defaulting to the [`FragmentVariant.Residue`](@ref FragmentVariant) variant.

# Supported keyword arguments
See [`Fragment`](@ref)
"""
@inline function Residue(ss::SecondaryStructure, number::Int; kwargs...)
    Fragment(ss, number; variant = FragmentVariant.Residue, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Returns `true` if the given fragment is a [`FragmentVariant.Residue`](@ref FragmentVariant),
`false` otherwise.
"""
@inline function isresidue(frag::Fragment)
    frag.variant === FragmentVariant.Residue
end

"""
    residue_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Fragment{T}

Returns the `Fragment{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
fragment exists or if the fragment is not a [`FragmentVariant.Residue`](@ref FragmentVariant).
"""
@inline function residue_by_idx(sys::System, idx::Int)
    frag = fragment_by_idx(sys, idx)
    isresidue(frag) || throw(KeyError(idx))
    frag
end

@inline function residue_by_idx(idx::Int)
    residue_by_idx(default_system(), idx)
end

"""
    residues(::Chain)
    residues(::ChainTable)
    residues(::SecondaryStructure)
    residues(::SecondaryStructureTable)
    residues(::Molecule)
    residues(::MoleculeTable)
    residues(::System = default_system())

Returns a `FragmentTable{T}` containing all [`FragmentVariant.Residue`](@ref FragmentVariant)
fragments of the given atom container or table.

# Supported keyword arguments
See [`fragments`](@ref)
"""
@inline function residues(
    ac::Union{Chain, ChainTable, SecondaryStructure, SecondaryStructureTable, Molecule, MoleculeTable, System} = default_system();
    kwargs...
)
    fragments(ac; variant = FragmentVariant.Residue, kwargs...)
end

"""
    nresidues(::Chain)
    nresidues(::SecondaryStructure)
    nresidues(::Molecule)
    nresidues(::System = default_system())

Returns the number of [`FragmentVariant.Residue`](@ref FragmentVariant) fragments in the given
atom container.

# Supported keyword arguments
See [`fragments`](@ref)
"""
@inline function nresidues(
    ac::Union{Chain, ChainTable, SecondaryStructure, SecondaryStructureTable, Molecule, MoleculeTable, System} = default_system();
    kwargs...
)
    nfragments(ac; variant = FragmentVariant.Residue, kwargs...)
end

"""
    nresidues(::ChainTable)
    nresidues(::SecondaryStructureTable)
    nresidues(::FragmentTable)
    nresidues(::MoleculeTable)

Returns the number of [`FragmentVariant.Residue`](@ref FragmentVariant) fragments in the given
table.
"""
@inline function nresidues(ft::FragmentTable)
    length(filter(isresidue, ft))
end

@inline function nresidues(ct::Union{ChainTable, SecondaryStructureTable, MoleculeTable})
    nfragments(ct; variant = FragmentVariant.Residue)
end

"""
    parent_residue(::Atom)

Returns the `Fragment{T}` containing the given atom. Returns `nothing` if no such fragment exists or
if the fragment is not a [`FragmentVariant.Residue`](@ref FragmentVariant).
"""
@inline function parent_residue(atom::Atom)
    frag = parent_fragment(atom)
    isnothing(frag) || !isresidue(frag) ? nothing : frag
end

function _fragment_variant(name::AbstractString)
    is_amino_acid(name) && return FragmentVariant.Residue
    is_nucleotide(name) && return FragmentVariant.Nucleotide
    FragmentVariant.None
end

# TODO adapt to variants
@inline function is_amino_acid(frag::Fragment)
    is_amino_acid(frag.name)
end

@inline function three_letter_code(frag::Fragment)
    @assert is_amino_acid(frag)
    frag.name
end

@inline function one_letter_code(frag::Fragment)
    @assert is_amino_acid(frag)
    one_letter_code(frag.name)
end

# TODO adapt to variants
@inline function is_nucleotide(frag::Fragment)
    is_nucleotide(frag.name)
end

# TODO: we should come up with a better test than just checking the name
is_nucleotide(name::AbstractString) = name ∈ ["A", "C", "G", "T", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI"]

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
