export AbstractAtomBijection, TrivialAtomBijection

abstract type AbstractAtomBijection{T} end

struct TrivialAtomBijection{T<:Real} <: AbstractAtomBijection{T}
    atoms_A::AbstractDataFrame
    atoms_B::AbstractDataFrame
 
    function TrivialAtomBijection{T}(A::AbstractMolecule{T}, B::AbstractMolecule{T}) where {T<:Real}
        new(atoms_df(A), atoms_df(B))
    end

    function TrivialAtomBijection{T}(atoms_A, B::AbstractMolecule{T}) where {T<:Real}
        atoms_B = atoms_df(B)
        new(atoms_A, atoms_B[atoms_B.number .∈ Ref(atoms_A.number), :])
    end
end

TrivialAtomBijection(A::AbstractMolecule, B::AbstractMolecule) = TrivialAtomBijection{Float32}(A, B)
TrivialAtomBijection(atoms_A, B::AbstractMolecule) = TrivialAtomBijection{Float32}(atoms_A, B)
