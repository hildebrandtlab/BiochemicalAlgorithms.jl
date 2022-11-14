using DataFrames

export PDBMolecule

mutable struct PDBMolecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    chains::Vector{PDBChain}

    properties::Properties

    function PDBMolecule{T}(name = "",
                        atoms = DataFrame(PDBAtom{T}[]),
                        bonds = DataFrame(Bond[]),
                        chains = PDBChain[],
                        properties = Properties()) where {T<:Real}
        new(name, atoms, bonds, chains, properties)
    end
end

PDBMolecule(name = "",
        atoms = DataFrame(PDBAtom{Float32}[]),
        bonds = DataFrame(Bond[]),
        chains = PDBChain[],
        properties = Properties()) = PDBMolecule{Float32}(name, atoms, bonds, chains, properties)
