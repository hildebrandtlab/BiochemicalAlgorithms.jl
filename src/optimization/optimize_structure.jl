export
    optimize_hydrogen_positions!,
    optimize_structure!

"""
    optimize_structure!(ff::ForceField)

Attempts to solve the energy optimization problem represented by the given force field object.

Internally `ff` is compiled into a cached, struct-of-arrays evaluator (see
[`compile`](@ref)) so the nonbonded pair list is reused across evaluations and
rebuilt only when atoms drift past a Verlet skin — instead of being rebuilt on
every energy/force evaluation. The canonical atom table (`atoms(system).r`/`.F`)
is kept current at every pair-list rebuild and fully synchronized when the
optimization finishes, so the rest of the library and any registered observables
see consistent state.

# Supported keyword arguments
 - `alg = OptimizationLBFGSB.LBFGSB()`
 - `accumulation::Type = T` — accumulation precision (`Float32`/`Float64`)
 - `backend::Symbol = :serial` — evaluation backend
 - `update_frequency::Integer = 0` — force a rebuild every N evaluations (0 = skin-only)
 - `max_displacement::Real = -1` — Verlet skin in Å (-1 = auto-derive from cutoffs)

All other keyword arguments are forwarded to
[Optimization.solve](https://docs.sciml.ai/Optimization/stable/API/solve/).
"""
function optimize_structure!(ff::ForceField{T};
        alg = OptimizationLBFGSB.LBFGSB(),
        accumulation::Type = T,
        backend::Symbol = :serial,
        update_frequency::Integer = 0,
        max_displacement::Real = -1,
        callback = nothing,
        kwargs...) where T
    
    cff = compile(ff; acc=accumulation, backend=backend,
                  update_frequency=update_frequency, max_displacement=max_displacement)

    _optimize_structure!(cff; alg=alg, callback=callback, kwargs...)
end

function optimize_structure!(cff::CompiledForceField{T};
        alg = OptimizationLBFGSB.LBFGSB(),
        callback = nothing,
        kwargs...) where T

    _optimize_structure!(cff; alg=alg, callback=callback, kwargs...)
end

function _optimize_structure!(cff::CompiledForceField{T};
        alg = OptimizationLBFGSB.LBFGSB(),
        callback = nothing,
        kwargs...) where T

    r0 = collect(Float64, Iterators.flatten(atoms(cff.ff.system).r))

    optf = Optimization.OptimizationFunction(
        (r, _ = nothing) -> begin
            prepare_eval!(cff, r)
            Float64(compute_energy!(cff))
        end,
        grad = (grad, r, _) -> begin
            prepare_eval!(cff, r)
            compute_forces!(cff)
            gradient_flat!(grad, cff)
            nothing
        end
    )
    prob = Optimization.OptimizationProblem(optf, r0)

    # Wrap any user callback so the canonical atom table is synchronized before
    # the callback fires (the Observable path notifies inside its callback, so
    # observers/visualization must see current coordinates and forces).
    solve_cb = callback === nothing ? nothing : (state, loss) -> begin
        sync_to_table!(cff; forces=true)
        callback(state, loss)
    end

    solution = solve_cb === nothing ?
        Optimization.solve(prob, alg; kwargs...) :
        Optimization.solve(prob, alg; callback=solve_cb, kwargs...)

    # write the final optimized state back into the canonical atom table
    set_positions_flat!(cff, solution.u)
    rebuild_pairlist!(cff)
    compute_forces!(cff)
    sync_to_table!(cff; forces=true)

    solution
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
