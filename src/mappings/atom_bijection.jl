export AbstractAtomBijection, TrivialAtomBijection

abstract type AbstractAtomBijection{T} end

struct TrivialAtomBijection{T<:Real} <: AbstractAtomBijection{T}
    atoms_A::DataFrame
    atoms_B::DataFrame
 
    function TrivialAtomBijection{T}(A::AbstractMolecule{T}, B::AbstractMolecule{T}) where {T<:Real}
        new(atoms(A), atoms(B))
    end

    function TrivialAtomBijection{T}(atoms_A, B::AbstractMolecule{T}) where {T<:Real}
        atoms_B = atoms(B)
        new(atoms_A, atoms_B[atoms_B.number .∈ Ref(atoms_A.number), :])
    end
end

TrivialAtomBijection(A::AbstractMolecule, B::AbstractMolecule) = TrivialAtomBijection{Float32}(A, B)
TrivialAtomBijection(atoms_A, B::AbstractMolecule) = TrivialAtomBijection{Float32}(atoms_A, B)