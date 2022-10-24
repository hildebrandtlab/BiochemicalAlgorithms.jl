using DataFrames

export Protein

mutable struct Protein{T<:Real} <: AbstractMolecule{T}
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
        atoms=DataFrame(PDBAtom[]), 
        bonds=DataFrame(Bond[]), 
        chains=DataFrame(ProteinChain[])) = Protein{Float32}(name, atoms, bonds, chains)
