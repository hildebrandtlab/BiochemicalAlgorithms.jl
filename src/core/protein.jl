using DataFrames

export Protein

mutable struct Protein{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    atom_properties::Dict{Int, AtomProperties} # mapping of atom.number => properties

    residues::DataFrame
    chains::DataFrame

    function Protein{T}(name = "", 
                        atoms = DataFrame(ProteinAtom{T}[]), 
                        bonds = DataFrame(Bond[]),
                        atom_properties=Dict{Int, AtomProperties}(),
                        residues = DataFrame(Residue[]),
                        chains = DataFrame(Chain[])) where {T<:Real}
        new(name, atoms, bonds, atom_properties, residues, chains)
    end
end

Protein(name="", 
        atoms=DataFrame(ProteinAtom[]), 
        bonds=DataFrame(Bond[]), 
        atom_properties=Dict{Int, AtomProperties}(),
        residues=DataFrame(Residue[]),
        chains=DataFrame(Chain[])) = Protein{Float32}(name, atoms, bonds, atom_properties, residues, chains)
