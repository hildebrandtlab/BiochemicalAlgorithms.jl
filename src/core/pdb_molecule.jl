using DataFrames

export PDBMolecule

mutable struct PDBMolecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    atom_properties::Dict{Int, AtomProperties} # mapping of atom.number => properties

    chains::Vector{PDBChain}

    function PDBMolecule{T}(name = "", 
                        atoms = DataFrame(PDBAtom{T}[]), 
                        bonds = DataFrame(Bond[]),
                        atom_properties=Dict{Int, AtomProperties}(),
                        chains = PDBChain[]) where {T<:Real}
        new(name, atoms, bonds, atom_properties, chains)
    end
end

PDBMolecule(name="", 
        atoms=DataFrame(PDBAtom[]), 
        bonds=DataFrame(Bond[]), 
        atom_properties=Dict{Int, AtomProperties}(),
        chains=PDBChain[]) = PDBMolecule{Float32}(name, atoms, bonds, atom_properties, chains)
