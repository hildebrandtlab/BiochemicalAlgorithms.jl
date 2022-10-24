using DataFrames

export AbstractMolecule, Molecule, count_atoms, count_bonds

abstract type AbstractMolecule{T} end

mutable struct Molecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame

    function Molecule{T}(name = "", 
                         atoms = DataFrame(Atom{T}[]), 
                         bonds = DataFrame(Bond[])) where {T<:Real}
        new(name, atoms, bonds)
    end
end

Molecule(name="", 
         atoms=DataFrame(Atom[]), 
         bonds=DataFrame(Bond[])) = Molecule{Float32}(name, atoms, bonds)

### Functions

function Base.push!(m::Molecule, atom::Atom)
    push!(m.atoms, atom)
end

function Base.push!(m::AbstractMolecule, bond::Bond)
    push!(m.bonds, bond)
end

function count_atoms(m::AbstractMolecule)
    # this is slightly complicated, because we need to unify by frame id
    # note that we assume that all frames have the same number of atoms

    # TODO: find a better way to handle this
    isempty(m.atoms) ? 0 : nrow(groupby(m.atoms, :frame_id)[1])
end

function count_bonds(m::AbstractMolecule)
    nrow(m.bonds)
end
