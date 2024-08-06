# removed in v0.3
@deprecate RigidTransform{T}(r::RotMatrix3{T}, t::Vector3{T}) where {T <: Real} RigidTransform(r, t)
@deprecate RigidTransform{T}(r::Matrix3{T}, t::Vector3{T}) where {T <: Real} RigidTransform(r, t)
@deprecate TrivialAtomBijection{T}(A::AbstractAtomContainer{T}, B::AbstractAtomContainer{T}) where T TrivialAtomBijection(A, B)
@deprecate TrivialAtomBijection{T}(A, B::AbstractAtomContainer{T}) where T TrivialAtomBijection(A, B)
@deprecate compute_energy(c::AbstractForceFieldComponent) compute_energy!(c)
@deprecate compute_energy(ff::ForceField) compute_energy!(ff)
@deprecate compute_forces compute_forces!
@deprecate compute_rmsd_minimizer(f::AbstractAtomBijection{T}, mini::Type{<: AbstractRMSDMinimizer}) where T compute_rmsd_minimizer(f; minimizer=mini)
