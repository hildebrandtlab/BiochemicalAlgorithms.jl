using DataFrames

export Protein

mutable struct Protein{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    residues::DataFrame
    chains::DataFrame

    function Protein{T}(name = "", 
                        atoms = DataFrame(ProteinAtom{T}[]), 
                        bonds = DataFrame(Bond[]),
                        residues = DataFrame(Residue[]),
                        chains = DataFrame(Chain[])) where {T<:Real}
        new(name, atoms, bonds, residues, chains)
    end
end

Protein(name="", 
        atoms=DataFrame(ProteinAtom[]), 
        bonds=DataFrame(Bond[]), 
        residues=DataFrame(Residue[])) = Protein{Float32}(name, atoms, bonds, residues)
