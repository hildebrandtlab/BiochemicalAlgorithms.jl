using AutoHashEquals
using DataFrames

export AbstractMolecule, Molecule, count_atoms, count_bonds, 
    has_property, get_property, set_property

abstract type AbstractMolecule{T} end

@auto_hash_equals mutable struct Molecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame
    properties::Properties

    function Molecule{T}(name = "",
                         atoms = DataFrame(Atom{T}[]),
                         bonds = DataFrame(Bond[]),
                         properties = Properties()) where {T<:Real}
        new(name, atoms, bonds, properties)
    end
end

Molecule(name = "",
         atoms = DataFrame(Atom{Float32}[]),
         bonds = DataFrame(Bond[]), 
         properties = Properties()) = Molecule{Float32}(name, atoms, bonds, properties)

### Functions

function Base.push!(m::Molecule, atom::Atom)
    push!(m.atoms, atom)
end

function Base.push!(m::AbstractMolecule, bond::Bond)
    push!(m.bonds, bond)
end

function count_atoms(m::AbstractMolecule)
    # this is slightly complicated, because we need to unify by frame id
    # note that we assume that all frames have the same number of atoms

    # TODO: find a better way to handle this
    isempty(m.atoms) ? 0 : nrow(groupby(m.atoms, :frame_id)[1])
end

function count_bonds(m::AbstractMolecule)
    nrow(m.bonds)
end

function has_property(m::AbstractMolecule, key::String)
    haskey(m.properties, key)
end

function get_property(m::AbstractMolecule, key::String)
    m.properties[key]
end

function set_property(m::AbstractMolecule, key::String, value)
    m.properties[key] = value
end