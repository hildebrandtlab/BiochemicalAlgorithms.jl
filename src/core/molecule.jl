using DataFrames

export AbstractMolecule, Molecule

abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

mutable struct Molecule{T<:Real} <: AbstractMolecule{T}
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
