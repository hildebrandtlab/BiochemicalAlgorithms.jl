export
    optimize_hydrogen_positions!,
    optimize_structure!

"""
    optimize_structure!(ff::ForceField)

Attempts to solve the energy optimization problem represented by the given force field object.

# Supported keyword arguments
This function passes all keyword arguments to
[Optimization.solve](https://docs.sciml.ai/Optimization/stable/API/solve/),
with the following default values:
 - `alg = Optimization.LBFGS()`
"""
function optimize_structure!(ff::ForceField; alg = Optimization.LBFGS(), kwargs...)
    r0 = collect(Float64, Iterators.flatten(atoms(ff.system).r))

    optf = Optimization.OptimizationFunction(
        (r, _ = nothing) -> begin
            atoms(ff.system).r .= eachcol(reshape(r, 3, :))
            update!(ff)
            compute_energy!(ff)
        end,
        grad = (grad, r, _) -> begin
            update!(ff)
            compute_forces!(ff)
            F = atoms(ff.system).F
            F[ff.constrained_atoms] .= Ref(zeros(3))
            grad .= -collect(Float64, Iterators.flatten(F)) ./ force_prefactor
            nothing
        end
    )
    prob = Optimization.OptimizationProblem(optf, r0)

    Optimization.solve(prob, alg; kwargs...)
end

"""
    optimize_hydrogen_positions!(ff::ForceField)

Variant of [`optimize_structure!`](@ref) that only optimizes hydrogen atom positions.

# Supported keyword arguments
Same as [`optimize_structure!`](@ref)
"""
function optimize_hydrogen_positions!(ff::ForceField; kwargs...)
    old_constrained = copy(ff.constrained_atoms)
    empty!(ff.constrained_atoms)
    append!(ff.constrained_atoms, findall(a -> a.element != Elements.H, atoms(ff.system)))
    solution = optimize_structure!(ff; kwargs...)
    empty!(ff.constrained_atoms)
    append!(ff.constrained_atoms, old_constrained)
    solution
end

for fun in [:optimize_structure!, :optimize_hydrogen_positions!]
    @eval begin
        function $(fun)(ff::Observables.Observable{ForceField{T}}; notification_frequency::Int = 1, kwargs...) where T
            iterations = 0
            _update = () -> begin
                iterations += 1
                iterations % notification_frequency == 0 && notify(ff)
                false
            end
            $(fun)(ff[]; callback=(θ, l) -> _update(), kwargs...)
        end
    end
end
