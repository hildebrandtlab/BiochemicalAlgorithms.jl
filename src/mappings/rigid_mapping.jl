using Statistics: mean
using LinearAlgebra: Hermitian, eigen

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
    DataFramesMeta.@with atoms_df(m) begin
       :r .= Ref(t) .+ :r
    end
    m
end

function rigid_transform!(m::AbstractMolecule{T}, transform::RigidTransform{T}) where {T<:Real}
    DataFramesMeta.@with atoms_df(m) begin
        :r .= Ref(transform.rotation) .* :r .+ Ref(transform.translation)
    end
    m
end

function compute_rmsd(f::AbstractAtomBijection{T}) where {T<:Real}
    r_BA = Vector{Vector3{T}}(undef, size(f.atoms_A, 1))
    DataFramesMeta.@with f.atoms_A r_BA .= :r
    DataFramesMeta.@with f.atoms_B r_BA .-= :r
    sqrt(mean(map(r -> transpose(r) * r, r_BA)))
end

function compute_rmsd_minimizer(f::AbstractAtomBijection{T}) where {T<:Real}
    r_A = Vector{Vector3{T}}(undef, size(f.atoms_A, 1))
    DataFramesMeta.@with f.atoms_A r_A .= :r
    mean_A = mean(r_A)

    r_B = Vector{Vector3{T}}(undef, size(f.atoms_B, 1))
    DataFramesMeta.@with f.atoms_B r_B .= :r
    mean_B = mean(r_B)

    R = mapreduce(t -> t[1] * transpose(t[2]), +, zip(r_B .- Ref(mean_B), r_A .- Ref(mean_A)))

    C = Hermitian(transpose(R) * R)
    μ, a = eigen(C)

    RigidTransform{T}(mapreduce(i -> 1/√μ[i] * (R * a[:, i]) * transpose(a[:, i]), +, 1:3), mean_B - mean_A)
end

compute_rmsd_minimizer(f) = compute_rmsd_minimizer{Float32}(f)

function map_rigid!(A::AbstractMolecule{T}, B::AbstractMolecule{T}; heavy_atoms_only::Bool = false) where {T<:Real}
    atoms_A = atoms_df(A)
    heavy_atoms_only && (atoms_A = atoms_A[atoms_A.element .!== Elements.H, :])

    U = compute_rmsd_minimizer(TrivialAtomBijection(atoms_A, B))

    rigid_transform!(A, U)

    A
end
