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
is not updated during the optimization, except right before a callback fires, and
is fully synchronized when the optimization finishes.

Force fields with a component `compile` does not support are optimized through
the slower reference path instead, which ignores `accumulation`, `backend`,
`update_frequency` and `max_displacement`.

After the optimization, `atoms(system).F` holds zero forces on constrained atoms.

# Callbacks
Before a callback fires, `atoms(system).r`/`.F` are set to the coordinates in
`state.u` and the forces at those coordinates, and `loss` is the total energy
reported by the optimizer for that state.

!!! warning
    Do **not** call `compute_energy!(ff)` or `compute_forces!(ff)` inside a
    callback (or any observer) to track energies during the optimization. The
    optimization evaluates the compiled copy of `ff`; the reference nonbonded
    caches of `ff` (stored pair distances) are not updated between iterations,
    so these calls would combine current bonded terms with stale nonbonded terms
    and return a wrong energy. Read `loss` instead. The reference
    force field is made consistent again once `optimize_structure!` returns.

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

    if !all(_is_compilable, ff.components)
        @warn "optimize_structure!: force field has components that cannot be compiled; " *
            "using the slower reference path (accumulation, backend, update_frequency " *
            "and max_displacement are ignored)." maxlog=1
        return _optimize_structure_reference!(ff, alg, callback; kwargs...)
    end

    cff = compile(ff; acc=accumulation, backend=backend,
                  update_frequency=update_frequency, max_displacement=max_displacement)

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

    # coordinates of the forces currently held in `cff.F`
    forces_at = fill(NaN, length(r0))

    optf = Optimization.OptimizationFunction(
        (r, _ = nothing) -> begin
            prepare_eval!(cff, r)
            Float64(compute_energy!(cff))
        end,
        grad = (grad, r, _) -> begin
            prepare_eval!(cff, r)
            compute_forces!(cff)
            copyto!(forces_at, r)
            gradient_flat!(grad, cff)
            nothing
        end
    )
    prob = Optimization.OptimizationProblem(optf, r0)

    # Wrap any user callback so the canonical atom table is synchronized before
    # the callback fires (the Observable path notifies inside its callback, so
    # observers/visualization must see current coordinates and forces). The last
    # evaluation need not have been at `state.u` (e.g. a rejected line-search
    # trial), so forces are recomputed unless they already belong to it.
    solve_cb = callback === nothing ? nothing : (state, loss) -> begin
        prepare_eval!(cff, state.u)
        if state.u != forces_at
            compute_forces!(cff)
            copyto!(forces_at, state.u)
        end
        sync_to_table!(cff; forces=true)
        callback(state, loss)
    end

    solution = solve_cb === nothing ?
        Optimization.solve(prob, alg; kwargs...) :
        Optimization.solve(prob, alg; callback=solve_cb, kwargs...)

    # write the final optimized state back into the canonical atom table
    set_positions_flat!(cff, solution.u)
    sync_to_table!(cff; forces=false)

    _refresh_reference_state!(ff)

    solution
end

# Leave the force field consistent with the optimized geometry for the reference
# path: `update!` refreshes each component's cached state (the nonbonded
# interaction list stores per-pair distances), and `compute_forces!` repopulates
# `atoms(system).F`. After this, the usual `compute_energy!(ff)` /
# `compute_forces!(ff)` report the minimized state.
function _refresh_reference_state!(ff::ForceField)
    update!(ff)
    compute_forces!(ff)
    _zero_constrained_forces!(ff)
    compute_energy!(ff)  # ensure the final energy is up-to-date
    nothing
end

function _zero_constrained_forces!(ff::ForceField{T}) where T
    atoms(ff.system).F[ff.constrained_atoms] .= Ref(zero(Vector3{T}))
    nothing
end

# Optimization through the components' own `compute_energy!`/`compute_forces!`,
# for force fields that `compile` cannot lower.
function _optimize_structure_reference!(ff::ForceField, alg, callback; kwargs...)
    r0 = collect(Float64, Iterators.flatten(atoms(ff.system).r))

    set_positions! = r -> begin
        atoms(ff.system).r .= eachcol(reshape(r, 3, :))
        update!(ff)
    end

    optf = Optimization.OptimizationFunction(
        (r, _ = nothing) -> begin
            set_positions!(r)
            Float64(compute_energy!(ff))
        end,
        grad = (grad, r, _) -> begin
            set_positions!(r)
            compute_forces!(ff)
            _zero_constrained_forces!(ff)
            grad .= -collect(Float64, Iterators.flatten(atoms(ff.system).F))
            nothing
        end
    )
    prob = Optimization.OptimizationProblem(optf, r0)

    solve_cb = callback === nothing ? nothing : (state, loss) -> begin
        set_positions!(state.u)
        compute_forces!(ff)
        _zero_constrained_forces!(ff)
        callback(state, loss)
    end

    solution = solve_cb === nothing ?
        Optimization.solve(prob, alg; kwargs...) :
        Optimization.solve(prob, alg; callback=solve_cb, kwargs...)

    atoms(ff.system).r .= eachcol(reshape(solution.u, 3, :))
    _refresh_reference_state!(ff)

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
    try
        optimize_structure!(ff; kwargs...)
    finally
        empty!(ff.constrained_atoms)
        append!(ff.constrained_atoms, old_constrained)
    end
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
