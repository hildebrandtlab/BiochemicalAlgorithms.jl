export
    AbstractRMSDMinimizer,
    RMSDMinimizerCoutsias,
    RMSDMinimizerKabsch,
    RigidTransform,
    compute_rmsd,
    compute_rmsd_minimizer,
    map_rigid!,
    rigid_transform!,
    translate!

"""
    $(TYPEDEF)

Abstract base type for RMSD minimizers.
"""
abstract type AbstractRMSDMinimizer end

"""
    $(TYPEDEF)

Abstract base type for RMSD minimizers based on <https://doi.org/10.1107/S0567739476001873>.
"""
abstract type RMSDMinimizerKabsch <: AbstractRMSDMinimizer end

"""
    $(TYPEDEF)

Abstract base type for RMSD minimizers based on <https://doi.org/10.1002/jcc.20110>.
"""
abstract type RMSDMinimizerCoutsias <: AbstractRMSDMinimizer end

"""
    $(TYPEDEF)

Rigid transformation represented by a single rotation (around a specific center of
rotation) and translation.

# Constructors
```julia
RigidTransform(r::RotMatrix3{T}, t::Vector3{T}, center::Vector3{T} = zeros(Vector3{T}))
RigidTransform(r::Matrix3{T}, t::Vector3{T}, center::Vector3{T} = zeros(Vector3{T}))
```
Creates a new `RigidTransform{T}` from the given rotation `r`, `center` of rotation, and
translation `t`.

!!! note
    From the documentation of Rotations.jl:
    > The given `Matrix3{T}` should have the property `I =RR^T`, but this isn't enforced
    > by the constructor.
"""
struct RigidTransform{T<:Real}
    rotation::RotMatrix3{T}
    translation::Vector3{T}
    center::Vector3{T}

    @inline function RigidTransform(
        r::RotMatrix3{T},
        t::Vector3{T},
        center::Vector3{T} = zeros(Vector3{T})
    ) where T
        new{T}(r, t, center)
    end

    @inline function RigidTransform(
        r::Matrix3{T},
        t::Vector3{T},
        center::Vector3{T} = zeros(Vector3{T})
    ) where T
        new{T}(RotMatrix3(r), t, center)
    end
end

"""
    translate!(::AtomTable{T}, t::Vector3{T})
    translate!(::AbstractAtomContainer{T}, t::Vector3{T})

Translates all atoms of the given container according to the given translation vector `t`.
"""
@inline function translate!(at::AtomTable{T}, t::Vector3{T}) where T
    at.r .+= Ref(t)
    at
end

@inline function translate!(ac::AbstractAtomContainer{T}, t::Vector3{T}) where T
    translate!(atoms(ac), t)
    ac
end

"""
    rigid_transform!(at::AtomTable{T}, transform::RigidTransform{T})
    rigid_transform!(ac::AbstractAtomContainer, transform::RigidTransform)

Applies the rotation and the translation represented by `transform` (in this order) to all
atoms of the given container.
"""
@inline function rigid_transform!(at::AtomTable{T}, transform::RigidTransform{T}) where T
    translate!(at, -transform.center)
    at.r .= Ref(transform.rotation) .* at.r .+ Ref(transform.translation) .+ Ref(transform.center)
    at
end

@inline function rigid_transform!(ac::AbstractAtomContainer, transform::RigidTransform)
    rigid_transform!(atoms(ac), transform)
    ac
end

"""
    compute_rmsd(f::AbstractAtomBijection)
    compute_rmsd(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T})
    compute_rmsd(A::AtomTable{T}, B::AtomTable{T})

Computes the [root mean square deviation (RMSD)](https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions)
of the given atom bijection. Defaults to [`TrivialAtomBijection`](@ref) if arguments are atom containers or tables.
"""
@inline function compute_rmsd(A::AtomTable{T}, B::AtomTable{T}) where T
    sqrt(mean(map(r -> transpose(r) * r, A.r .- B.r)))
end

@inline function compute_rmsd(f::AbstractAtomBijection)
    compute_rmsd(atoms(f)...)
end

@inline function compute_rmsd(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}) where T
    compute_rmsd(TrivialAtomBijection(A, B))
end

"""
    compute_rmsd_minimizer(f::AbstractAtomBijection)
    compute_rmsd_minimizer(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T})
    compute_rmsd_miminizer(A::AtomTable{T}, B::AtomTable{T})

Computes the RMSD-minimizing rigid transformation for the given atom bijection. Defaults
to [`TrivialAtomBijection`](@ref) if arguments are atom containers or tables.

# Supported keyword arguments
 - `minimizer::Type{<: AbstractRMSDMinimizer} = RMSDMinimizerCoutsias`
   Method used for computing the optimal rotation matrix. See [`RMSDMinimizerCoutsias`](@ref)
   and [`RMSDMinimizerKabsch`](@ref).
"""
function compute_rmsd_minimizer(
    A::AtomTable{T},
    B::AtomTable{T};
    minimizer::Type{<:AbstractRMSDMinimizer} = RMSDMinimizerCoutsias
) where T
    mean_A = mean(A.r)
    mean_B = mean(B.r)

    R = mapreduce(t -> t[1] * transpose(t[2]), +, zip(A.r .- Ref(mean_A), B.r .- Ref(mean_B)))

    RigidTransform(_compute_rotation(R, minimizer), mean_B .- mean_A, mean_A)
end

@inline function compute_rmsd_minimizer(f::AbstractAtomBijection; kwargs...)
    compute_rmsd_minimizer(atoms(f)...; kwargs...)
end

@inline function compute_rmsd_minimizer(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}; kwargs...) where T
    compute_rmsd_minimizer(TrivialAtomBijection(A, B); kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Computes the rotation matrix by solving the eigenvalue problem given as the correlation matrix `C`.
Uses all resulting eigenvalues and eigenvectors.
Warns if the correlation matrix is not positive definit (contains negative eigenvalues or eigenvalues equal to 0)
and uses the alternative approch [`RMSDMinimizerCoutsias`](@ref) instead.
Returns a `RotMatrix3`.
"""
function _compute_rotation(R::Matrix3{T}, ::Type{RMSDMinimizerKabsch}) where {T<:Real}
    C = Hermitian(transpose(R) * R)
    μ, a = eigen(C)

    if minimum(μ) <= 0
        @warn("Correlation matrix is not positive definit! Computing rotation through `RMSDMinimizerCoutsias` instead!")
        return _compute_rotation(R, RMSDMinimizerCoutsias)
    end

    RotMatrix3{T}(mapreduce(i -> 1/√μ[i] * (R * a[:, i]) * transpose(a[:, i]), +, 1:3))
end

"""
    $(TYPEDSIGNATURES)

Computes the rotation matrix by solving the eigenvalue problem given as the residual matrix `F`.
Uses only the largest of the resulting eigenvalues to generate the quaternion describing the
optimal rotation that maps the atoms onto each other.
Returns a `RotMatrix3`.
"""
function _compute_rotation(R::Matrix3{T}, ::Type{RMSDMinimizerCoutsias}) where {T<:Real}
    # Residual matrix F
    F = zeros(4,4)
    F[1,1] = R[1,1] + R[2,2] + R[3,3]
    F[2,1] = R[2,3] - R[3,2]
    F[3,1] = R[3,1] - R[1,3]
    F[4,1] = R[1,2] - R[2,1]

    F[1,2] = R[2,3] - R[3,2]
    F[2,2] = R[1,1] - R[2,2] - R[3,3]
    F[3,2] = R[1,2] + R[2,1]
    F[4,2] = R[1,3] + R[3,1]

    F[1,3] = R[3,1] - R[1,3]
    F[2,3] = R[1,2] + R[2,1]
    F[3,3] = -R[1,1] + R[2,2] - R[3,3]
    F[4,3] = R[2,3] + R[3,2]

    F[1,4] = R[1,2] - R[2,1]
    F[2,4] = R[1,3] + R[3,1]
    F[3,4] = R[2,3] + R[3,2]
    F[4,4] = -R[1,1] - R[2,2] + R[3,3]

    μ, a = eigen(F)
    i = argmax(μ)

    RotMatrix3{T}(QuatRotation(quat(a[1,i], a[2,i], a[3,i], a[4,i])))
end

"""
    map_rigid!(f::AbstractAtomBijection)
    map_rigid!(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T})
    map_rigid!(A::AtomTable{T}, B::AtomTable{T})

Computes, applies, and returns the RMSD-minimizing rigid transformation for the given atom bijection.
Defaults to [`TrivialAtomBijection`](@ref) if arguments are atom containers or tables. Only the first
structure (`A`) is modified, i.e. mapped onto the second one (`B`).

# Supported keyword arguments
 - `heavy_atoms_only::Bool = false`
   If `true`, hydrogen atoms are ignored during the computation of the optimal transformation. Otherwise,
   all atoms are used.
 - `minimizer::Type{<: AbstractRMSDMinimizer} = RMSDMinimizerCoutsias`
   See [`compute_rmsd_minimizer`](@ref)
"""
function map_rigid!(
    A::AtomTable{T},
    B::AtomTable{T};
    heavy_atoms_only::Bool = false,
    minimizer::Type{<: AbstractRMSDMinimizer} = RMSDMinimizerCoutsias
) where T
    refA = heavy_atoms_only ? filter(atom -> atom.element != Elements.H, A) : A
    refB = heavy_atoms_only ? filter(atom -> atom.element != Elements.H, B) : B

    rt = compute_rmsd_minimizer(refA, refB; minimizer = minimizer)
    rigid_transform!(A, rt)
    rt
end

@inline function map_rigid!(f::AbstractAtomBijection; kwargs...)
    map_rigid!(atoms(f)...; kwargs...)
end

@inline function map_rigid!(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}; kwargs...) where T
    map_rigid!(TrivialAtomBijection(A, B); kwargs...)
end

"""
    $(TYPEDSIGNATURES)

The transformation maps
 1. the point `w1` onto the point `v1` and
 2. the point `w2` onto the ray that starts in `v1` and goes through `v2`
 3. the point `w3` into the plane generated by `v1`, `v2` and `v3`
"""
function match_points(
        w1::Vector3{T}, w2::Vector3{T}, w3::Vector3{T},
        v1::Vector3{T}, v2::Vector3{T}, v3::Vector3{T}) where {T<:Real}
    ϵ = T(0.00001)
    ϵ₂ = T(0.00000001)

    # Compute the translations that map v1 and w1 onto the origin 
    # and apply them to v2, v3 and w2, w3.
    tw2 = w2 - w1
    tw3 = w3 - w1

    tv2 = v2 - v1
    tv3 = v3 - v1

    dist_v2_v1 = squared_norm(tv2)
    dist_w2_w1 = squared_norm(tw2)
    dist_w3_w1 = squared_norm(tw3)
    dist_v3_v1 = squared_norm(tv3)

    # Try to remove nasty singularities arising if the first two
    # points in each point set are too close to each other:
    #   (a) ensure (v2 != v1) 
    if ((dist_v2_v1 < ϵ₂) && (dist_v3_v1 >= ϵ₂))
        tv3, tv2 = tv2, tv3
    end

    #   (b) ensure (w2 != w1) 
    if ((dist_w2_w1 < ϵ₂) && (dist_w3_w1 >= ϵ₂))
        tw3, tw2 = tw2, tw3
    end

    # initialize translation
    final_translation = -w1
    final_rotation = T(1)I(3)

    if ((squared_norm(tv2) >= ϵ₂) && (squared_norm(tw2) >= ϵ₂))
        # calculate the rotation axis: orthogonal to tv2 and tw2
        tw2 = normalize(tw2)
        tv2 = normalize(tv2)

        rotation_axis = tw2 + tv2

        rotation = if (squared_norm(rotation_axis) < ϵ)
            # the two axes seem to be antiparallel -
            # invert the second vector
            T(-1)I(3)
        else
            # rotate around the rotation axis
            AngleAxis{T}(π, rotation_axis...)
        end

        tw2 = rotation * tw2
        tw3 = rotation * tw3

        final_rotation    = rotation * final_rotation
        final_translation = rotation * final_translation

        if ((squared_norm(tw3) > ϵ₂) && (squared_norm(tv3) > ϵ₂))
            tw3 = normalize(tw3)
            tv3 = normalize(tv3)

            axis_w = cross(tv2, tw3)
            axis_v = cross(tv2, tv3)

            if ((squared_norm(axis_v) > ϵ₂) && (squared_norm(axis_w) > ϵ₂))
                axis_v = normalize(axis_v)
                axis_w = normalize(axis_w)

                rotation_axis = cross(axis_w, axis_v)

                if (squared_norm(rotation_axis) < ϵ₂)
                    scalar_prod = dot(axis_w, axis_v)
                    rotation = if (scalar_prod < 0.0)
                        AngleAxis{T}(π, tv2...)
                    else
                        T(1)I(3)
                    end
                else
                    # Compute the rotation that maps tw3 onto tv3
                    product = dot(axis_w, axis_v)
                    product = min(T(1.0), max(T(-1.0), product))

                    angle = acos(product)
                    rotation = if (angle > ϵ)
                        AngleAxis{T}(angle, rotation_axis...)
                    else
                        # Use the identity matrix instead.
                        T(1.0)I(3)
                    end
                end

                final_rotation    = rotation * final_rotation
                final_translation = rotation * final_translation
            end
        end
    end

    # apply the translation onto v1
    final_translation += v1

    # done
    return final_translation, final_rotation
end
