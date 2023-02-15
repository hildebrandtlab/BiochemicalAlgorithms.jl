export AbstractAtomContainer, atoms, bonds, chains, count_atoms, count_bonds, count_fragments, 
    count_nucleotides, count_residues, fragments, name, nucleotides, properties, residues, set_name!

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T <: Real} end

"""
    $(TYPEDSIGNATURES)

Returns a `DataFrame` containing all atoms in the given container.
"""
@inline function atoms(ac::C) where {T, C <: AbstractAtomContainer{T}}
    ac.atoms::DataFrame
end

"""
    $(TYPEDSIGNATURES)

Returns a `DataFrame` containing all bonds in the given container.
"""
@inline function bonds(ac::C) where {T, C <: AbstractAtomContainer{T}}
    ac.bonds::DataFrame
end

"""
    $(TYPEDSIGNATURES)

Returns a `DataFrame` containing all fragments in the given container.
"""
@inline function fragments(::C) where {T, C <: AbstractAtomContainer{T}}
    DataFrame(Fragment[])
end

"""
    $(TYPEDSIGNATURES)

Returns the name of the given container.
"""
@inline function name(ac::C) where {T, C <: AbstractAtomContainer{T}}
    ac.name::String
end

"""
    $(TYPEDSIGNATURES)

Returns a `DataFrame` containing all nucleotides in the given container.
"""
@inline function nucleotides(::C) where {T, C <: AbstractAtomContainer{T}}
    DataFrame(Nucleotide[])
end

"""
    $(TYPEDSIGNATURES)

Returns the properties dictionary of the given container.
"""
@inline function properties(ac::C) where {T, C <: AbstractAtomContainer{T}}
    ac.properties::Properties
end

"""
    $(TYPEDSIGNATURES)

Returns a `DataFrame` containing all residues in the given container.
"""
@inline function residues(::C) where {T, C <: AbstractAtomContainer{T}}
    DataFrame(Residue[])
end

"""
    $(TYPEDSIGNATURES)

Renames the given container.
"""
@inline function set_name!(ac::C, name::String) where {T, C <: AbstractAtomContainer{T}}
    ac.name::String = name
end

#=
    TODO document
=#
function Base.push!(ac::C, atom::Atom{T}) where {T, C <: AbstractAtomContainer{T}}
    push!(atoms(ac), atom)
end

#=
    TODO document
=#
function Base.push!(ac::C, bond::Bond) where {T, C <: AbstractAtomContainer{T}}
    push!(bonds(ac), bond)
end

"""
    $(TYPEDSIGNATURES)

Returns the total number of atoms in the given container.
"""
function count_atoms(ac::C) where {T, C <: AbstractAtomContainer{T}}
    nrow(atoms(ac))
end

"""
    $(TYPEDSIGNATURES)

Returns the total number of bonds in the given container.
"""
@inline function count_bonds(ac::C) where {T, C <: AbstractAtomContainer{T}}
    nrow(bonds(ac))
end

"""
    $(TYPEDSIGNATURES)

Returns the total number of fragments in the given container.
"""
@inline function count_fragments(ac::C) where {T, C <: AbstractAtomContainer{T}}
    nrow(fragments(ac))
end

"""
    $(TYPEDSIGNATURES)

Returns the total number of nucleotides in the given container.
"""
@inline function count_nucleotides(ac::C) where {T, C <: AbstractAtomContainer{T}}
    nrow(nucleotides(ac))
end

"""
    $(TYPEDSIGNATURES)

Returns the total number of residues in the given container.
"""
@inline function count_residues(ac::C) where {T, C <: AbstractAtomContainer{T}}
    nrow(residues(ac))
end
