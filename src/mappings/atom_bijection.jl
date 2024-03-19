export AbstractAtomBijection, TrivialAtomBijection

abstract type AbstractAtomBijection{T} end

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

TrivialAtomBijection(A::AbstractAtomContainer, B::AbstractAtomContainer) = TrivialAtomBijection{Float32}(A, B)
TrivialAtomBijection(atoms_A, B::AbstractAtomContainer) = TrivialAtomBijection{Float32}(atoms_A, B)
