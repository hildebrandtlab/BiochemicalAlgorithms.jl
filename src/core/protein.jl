using AutoHashEquals
using DataFrames

export Protein

@auto_hash_equals mutable struct Protein{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    chains::Vector{ProteinChain}

    function Protein{T}(name = "", 
                        atoms = DataFrame(PDBAtom{T}[]), 
                        bonds = DataFrame(Bond[]),
                        chains = ProteinChain[]) where {T<:Real}
        new(name, atoms, bonds, chains)
    end
end

Protein(name="",
        atoms=DataFrame(PDBAtom{Float32}[]),
        bonds=DataFrame(Bond[]),
        chains=ProteinChain[]) = Protein{Float32}(name, atoms, bonds, chains)
