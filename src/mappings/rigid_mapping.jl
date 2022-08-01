using Statistics: mean
using LinearAlgebra: Hermitian, eigvals, eigvecs

export AbstractRMSDMinimizer, RigidTransform, rigid_transform!, compute_rmsd_minimizer, compute_rmsd, translate!, map_rigid!

abstract type AbstractRMSDMinimizer end
abstract type RMSDMinimizerKabsch <: AbstractRMSDMinimizer end

struct RigidTransform{T<:Real}
    rotation::Matrix3{T}
    translation::Vector3{T}

    function RigidTransform{T}(r::Matrix3{T}, t::Vector3{T}) where {T<:Real}
        new(r, t)
    end
end

RigidTransform(r::Matrix3, t::Vector3) = RigidTransform{Float32}(r, t)


### Functions

function translate!(m::AbstractMolecule{T}, t::Vector3{T}) where {T<:Real}
    m.atoms.r = Ref(t) .+ m.atoms.r

    m
end

function rigid_transform!(m::AbstractMolecule{T}, transform::RigidTransform{T}) where {T<:Real}
    m.atoms.r = Ref(transform.rotation) .* m.atoms.r .+ Ref(transform.translation)

    m
end

function compute_rmsd(f::AbstractAtomBijection{T}) where {T<:Real}
    sqrt(mean(map(r -> transpose(r) * r, f.atoms_A.r .- f.atoms_B.r)))
end

function compute_rmsd_minimizer(f::AbstractAtomBijection{T}) where {T<:Real}
    r_A = f.atoms_A.r
    r_B = f.atoms_B.r

    mean_A = mean(r_A)
    mean_B = mean(r_B)

    R = mapreduce(t -> t[1] * transpose(t[2]), +, zip(r_B .- Ref(mean_B), r_A .- Ref(mean_A)))

    C = Hermitian(transpose(R) * R)

    μ = eigvals(C)
    a = eigvecs(C)

    RigidTransform{T}(mapreduce(i -> 1/√μ[i] * (R * a[:, i]) * transpose(a[:, i]), +, 1:3), mean_B - mean_A)
end

compute_rmsd_minimizer(f) = compute_rmsd_minimizer{Float32}(f)

function map_rigid!(A::AbstractMolecule{T}, B::AbstractMolecule{T}; heavy_atoms_only::Bool = false) where {T<:Real}
    atoms_A = !heavy_atoms_only ? A.atoms : A.atoms[A.atoms.element .!== Elements.H, :]

    U = compute_rmsd_minimizer(TrivialAtomBijection(atoms_A, B))

    rigid_transform!(A, U)

    A
end