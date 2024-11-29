export
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
        grad = (F, r, _) -> begin
            update!(ff)
            compute_forces!(ff)
            F .= -collect(Float64, Iterators.flatten(atoms(ff.system).F)) ./ force_prefactor
            nothing
        end
    )
    prob = Optimization.OptimizationProblem(optf, r0)

    Optimization.solve(prob, alg; kwargs...)
end
