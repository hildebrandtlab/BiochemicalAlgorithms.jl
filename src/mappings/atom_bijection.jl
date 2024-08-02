export
    AbstractAtomBijection,
    TrivialAtomBijection

"""
    $(TYPEDEF)

Abstract base type for atom bijections.
"""
abstract type AbstractAtomBijection{T} end

"""
    $(TYPEDSIGNATURES)

Returns the tuple of atom tables represented by this bijection.
"""
function atoms(ab::AbstractAtomBijection)
    error("atoms() not implemented for atom bijection type $(typeof(ab))")
end

"""
    TrivialAtomBijection{T} <: AbstractAtomBijection{T}

Bijection of atoms based on their order in the associated atom containers.

# Constructors
```julia
TrivialAtomBijection(A::AtomTable{T}, B::AtomTable{T})
```
Creates a new `TrivialAtomBijection{T}` from the given tables. Atom order is determined by
row number.

```julia
TrivialAtomBijection(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T})
```
Creates a new `TrivialAtomBijection{T}` with the atoms of the individual atom containers.

```julia
TrivialAtomBijection(A::AtomTable{T}, B::AbstractAtomContainer{T})
```
Creates a new `TrivialAtomBijection{T}` by filtering `B` for atom numbers contained in `A`.
Atom order after filtering must be the same as in `A`.
"""
struct TrivialAtomBijection{T} <: AbstractAtomBijection{T}
    atoms_A::AtomTable{T}
    atoms_B::AtomTable{T}
 
    @inline function TrivialAtomBijection(A::AtomTable{T}, B::AtomTable{T}) where T
        new{T}(A, B)
    end

    @inline function TrivialAtomBijection(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}) where T
        new{T}(atoms(A), atoms(B))
    end

    @inline function TrivialAtomBijection(A::AtomTable{T}, B::AbstractAtomContainer{T}) where T
        anum = Set(A.number)
        new{T}(A, filter(atom -> atom.number in anum, atoms(B)))
    end
end

@inline function atoms(ab::TrivialAtomBijection)
    (ab.atoms_A, ab.atoms_B)
end
