export AbstractAtomBijection, TrivialAtomBijection

abstract type AbstractAtomBijection{T} end

struct TrivialAtomBijection{T} <: AbstractAtomBijection{T}
    atoms_A::AtomTable{T}
    atoms_B::AtomTable{T}
 
    function TrivialAtomBijection{T}(A::AbstractMolecule{T}, B::AbstractMolecule{T}) where T
        new(atoms(A), atoms(B))
    end

    function TrivialAtomBijection{T}(atoms_A, B::AbstractMolecule{T}) where T
        atoms_B = atoms(B)
        anum = Set(atoms_A.number)
        new(atoms_A, filter(atom -> atom.number in anum, atoms_B))
    end
end

TrivialAtomBijection(A::AbstractMolecule, B::AbstractMolecule) = TrivialAtomBijection{Float32}(A, B)
TrivialAtomBijection(atoms_A, B::AbstractMolecule) = TrivialAtomBijection{Float32}(atoms_A, B)
