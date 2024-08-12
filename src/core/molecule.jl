export
    Protein,
    Molecule,
    MoleculeTable,
    isprotein,
    molecule_by_idx,
    molecules,
    nmolecules,
    nproteins,
    parent_molecule,
    parent_protein,
    protein_by_idx,
    proteins

"""
    Molecule{T} <: AbstractAtomContainer{T}

Mutable representation of an individual molecule in a system.

# Public fields
 - `idx::Int`
 - `name::String`

# Private fields
 - `variant::MoleculeVariantType`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Molecule(
    sys::System{T};
    # keyword arguments
    name::String = "",
    variant::MoleculeVariantType = MoleculeVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{T}` in the given system.

```julia
Molecule(;
    #keyword arguments
    name::String = "",
    variant::MoleculeVariantType = MoleculeVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{Float32}` in the default system.
"""
const Molecule{T} = AtomContainer{T, _MoleculeTable}

@inline function Molecule(
    sys::System{T};
    kwargs...
) where T
    idx = _next_idx(sys)
    push!(sys._molecules, idx; kwargs...)
    molecule_by_idx(sys, idx)
end

@inline function Molecule(; kwargs...)
    Molecule(default_system(); kwargs...)
end

"""
    MoleculeTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system molecules (or a subset thereof). Molecule tables can be
generated using [`molecules`](@ref) or filtered from other molecule tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `name::AbstractVector{String}`

# Private columns
 - `variant::AbstractVector{MoleculeVariantType}`
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
"""
const MoleculeTable{T} = SystemComponentTable{T, Molecule{T}}

@inline function _filter_molecules(f::Function, sys::System{T}) where T
    MoleculeTable{T}(sys, _filter_idx(f, sys._molecules))
end

@inline _table(sys::System{T}, ::Type{Molecule{T}}) where T = sys._molecules

@inline function _hascolumn(::Type{<: Molecule}, nm::Symbol)
    nm in _molecule_table_cols_set || nm in _molecule_table_cols_priv
end

@doc raw"""
    parent_molecule(::Atom)
    parent_molecule(::Chain)
    parent_molecule(::Fragment)

Returns the `Molecule{T}` containing the given object. Returns `nothing` if no such molecule exists.
""" parent_molecule

"""
    push!(::System{T}, ::Molecule{T})

Creates a copy of the given molecule in the given system. The new molecule is automatically assigned a
new `idx`.
"""
@inline function Base.push!(sys::System{T}, mol::Molecule{T}) where T
    Molecule(sys;
        name = mol.name,
        variant = mol.variant,
        properties = mol.properties,
        flags = mol.flags
    )
    sys
end

"""
    molecule_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Molecule{T}

Returns the `Molecule{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
molecule exists.
"""
@inline function molecule_by_idx(sys::System{T}, idx::Int) where T
    Molecule{T}(sys, _row_by_idx(sys._molecules, idx))
end

@inline function molecule_by_idx(idx::Int)
    molecule_by_idx(default_system(), idx)
end

"""
    molecules(::System{T} = default_system())

Returns a `MoleculeTable{T}` containing all molecules of the given system.

# Supported keyword arguments
 - `variant::Union{Nothing, MoleculeVariantType} = nothing`
The keyword argument limits the results to molecules matching the given variant type.
Keyword arguments set to `nothing` are ignored.
"""
@inline function molecules(sys::System{T} = default_system();
    variant::Union{Nothing, MoleculeVariantType} = nothing
) where T
    isnothing(variant) && return MoleculeTable{T}(sys, sys._molecules.idx)
    _filter_molecules(mol -> mol.variant == something(variant), sys)
end

"""
    nmolecules(::System = default_system())

Returns the number of molecules in the given system.

# Supported keyword arguments
See [`molecules`](@ref)
"""
function nmolecules(sys::System = default_system(); kwargs...)
    length(molecules(sys; kwargs...))
end

#=
    Molecule atoms
=#
@inline atoms(mol::Molecule; kwargs...) = atoms(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline natoms(mol::Molecule; kwargs...) = natoms(parent(mol); molecule_idx = mol.idx, kwargs...)

@inline function Atom(mol::Molecule, number::Int, element::ElementType; kwargs...)
    Atom(parent(mol), number, element; molecule_idx = mol.idx, kwargs...)
end

@inline function Base.push!(mol::Molecule{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(mol), atom; molecule_idx = mol.idx, kwargs...)
    mol
end

#=
    Molecule bonds
=#
@inline bonds(mol::Molecule; kwargs...) = bonds(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline nbonds(mol::Molecule; kwargs...) = nbonds(parent(mol); molecule_idx = mol.idx, kwargs...)

#=
    Variant: Protein
=#
"""
    Protein(sys::System = default_system())

`Molecule` constructor defaulting to the [`MoleculeVariant.Protein`](@ref MoleculeVariant) variant.

# Supported keyword arguments
See [`Molecule`](@ref)
"""
@inline function Protein(sys::System = default_system(); kwargs...)
    Molecule(sys; variant = MoleculeVariant.Protein, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Returns `true` if the given molecule is a [`MoleculeVariant.Protein`](@ref MoleculeVariant),
`false` otherwise.
"""
@inline function isprotein(mol::Molecule)
    mol.variant === MoleculeVariant.Protein
end

"""
    protein_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Molecule{T}

Returns the `Molecule{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
molecule exists or if the molecule is not a [`MoleculeVariant.Protein`](@ref MoleculeVariant).
"""
@inline function protein_by_idx(sys::System, idx::Int)
    mol = molecule_by_idx(sys, idx)
    isprotein(mol) || throw(KeyError(idx))
    mol
end

@inline function protein_by_idx(idx::Int)
    protein_by_idx(default_system(), idx)
end

"""
    proteins(::System{T} = default_system())

Returns a `MoleculeTable{T}` containing all [`MoleculeVariant.Protein`](@ref MoleculeVariant) molecules
of the given system.

# Supported keyword arguments
See [`molecules`](@ref)
"""
@inline function proteins(sys::System = default_system(); kwargs...)
    molecules(sys; variant = MoleculeVariant.Protein, kwargs...)
end

"""
    nproteins(::System = default_system())

Returns the number of [`MoleculeVariant.Protein`](@ref MoleculeVariant) molecules in the given system.

# Supported keyword arguments
See [`molecules`](@ref)
"""
@inline function nproteins(sys::System = default_system(); kwargs...)
    nmolecules(sys; variant = MoleculeVariant.Protein, kwargs...)
end

"""
    parent_protein(::Atom)
    parent_protein(::Chain)
    parent_protein(::Fragment)

Returns the `Molecule{T}` containing the given object. Returns `nothing` if no such molecule exists or
if the molecule is not a [`MoleculeVariant.Protein`](@ref MoleculeVariant).
"""
@inline function parent_protein(ac)
    mol = parent_molecule(ac)
    isnothing(mol) || !isprotein(mol) ? nothing : mol
end
