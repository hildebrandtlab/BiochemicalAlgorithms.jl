using DataFrames

export Protein

mutable struct Protein{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    atom_properties::Dict{Int, AtomProperties} # mapping of atom.number => properties

    chains::Vector{ProteinChain}

    function Protein{T}(name = "", 
                        atoms = DataFrame(PDBAtom{T}[]), 
                        bonds = DataFrame(Bond[]),
                        atom_properties=Dict{Int, AtomProperties}(),
                        chains = ProteinChain[]) where {T<:Real}
        new(name, atoms, bonds, atom_properties, chains)
    end
end

Protein(name = "", 
        atoms = DataFrame(PDBAtom[]), 
        bonds = DataFrame(Bond[]), 
        atom_properties = Dict{Int, AtomProperties}(),
        chains = ProteinChain[]) = Protein{Float32}(name, atoms, bonds, atom_properties, chains)
