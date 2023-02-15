using DataFrames

export AbstractMolecule, Molecule

abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

mutable struct Molecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame
    properties::Properties

    function Molecule{T}(name = "") where T
        new(name, DataFrame(Atom{T}[]), DataFrame(Bond[]), Properties())
    end
end

Molecule(name = "") = Molecule{Float32}(name)
