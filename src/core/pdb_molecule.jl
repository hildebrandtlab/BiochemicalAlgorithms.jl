using DataFrames

export PDBMolecule

mutable struct PDBMolecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    chains::Vector{PDBChain}

    function PDBMolecule{T}(name = "",
                        atoms = DataFrame(PDBAtom{T}[]),
                        bonds = DataFrame(Bond[]),
                        chains = PDBChain[]) where {T<:Real}
        new(name, atoms, bonds, chains)
    end
end

PDBMolecule(name = "",
        atoms = DataFrame(PDBAtom{Float32}[]),
        bonds = DataFrame(Bond[]),
        chains = PDBChain[]) = PDBMolecule{Float32}(name, atoms, bonds, chains)
