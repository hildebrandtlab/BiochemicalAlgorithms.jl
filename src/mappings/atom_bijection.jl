export AbstractAtomBijection, TrivialAtomBijection
"""
    $(TYPEDEF)

Abstract base type for AtomBijections. 
    
"""
abstract type AbstractAtomBijection{T} end

"""
    TrivialAtomBijection{T} <: AbstractAtomBijection{T}
Mutable representation of a bijection of atoms based on the order of atoms in the individual atom containers, i.e., atom 1 of con

# Public fields
- `atoms_A::AtomTable{T}`
- `atoms_B::AtomTable{T}`

# Constructors
`TrivialAtomBijection{T}(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}) where T`

Creates a new `TrivialAtomBijection{T}` with the atoms of the individual atom containers. 

`TrivialAtomBijection{T}(atoms_A, B::AbstractAtomContainer{T}) where T`
Creates a new `TrivialAtomBijection{T}`based on the unique set of atom_numbers of `atoms_A`.
"""
struct TrivialAtomBijection{T} <: AbstractAtomBijection{T}
    atoms_A::AtomTable{T}
    atoms_B::AtomTable{T}
 
    function TrivialAtomBijection{T}(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}) where T
        new(atoms(A), atoms(B))
    end

    function TrivialAtomBijection{T}(atoms_A, B::AbstractAtomContainer{T}) where T
        atoms_B = atoms(B)
        anum = Set(atoms_A.number)
        new(atoms_A, filter(atom -> atom.number in anum, atoms_B))
    end
end

"""
    TrivialAtomBijection(A::AbstractAtomContainer, B::AbstractAtomContainer)

Creates a new `TrivialAtomBijection{Float32}` as a default.
"""
TrivialAtomBijection(A::AbstractAtomContainer, B::AbstractAtomContainer) = TrivialAtomBijection{Float32}(A, B)

"""
    TrivialAtomBijection(A::AbstractAtomContainer, B::AbstractAtomContainer)

Creates a new `TrivialAtomBijection{Float32}` based on the unique set of atom_numbers of `atoms_A`.
"""
TrivialAtomBijection(atoms_A, B::AbstractAtomContainer) = TrivialAtomBijection{Float32}(atoms_A, B)
