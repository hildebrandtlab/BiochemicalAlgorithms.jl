export optimize_structure!, optimize_hydrogen_positions!

function compute_energy(ff::ForceField{T}, r::AbstractVector{T}; verbose=false) where {T<:Real}
    # first, convert the vector of positions to something we can use
    atoms_df(ff.system).r = Vector3{T}.(r[1:3:end], r[2:3:end], r[3:3:end])

    update!(ff)

    # then, compute the energy
    compute_energy(ff; verbose)
end

function frule((_, _, Δr), ::typeof(compute_energy), ff, r)
    energy = compute_energy(ff, r)
    compute_forces(ff)
    grad = -(atoms_df(ff.system).F / force_prefactor)
    grad[ff.constrained_atoms] .= Ref(zero(Vector3{T}))
    return energy, grad * Δr
end

function rrule(::typeof(compute_energy), ff::ForceField{T}, r) where {T<:Real}
    function compute_energy_pullback(Δr)
        compute_forces(ff)
        grad = -(atoms_df(ff.system).F / force_prefactor)
        grad[ff.constrained_atoms] .= Ref(zero(Vector3{T}))

        (NoTangent(), NoTangent(), collect(Iterators.flatten(grad)) * Δr)
    end

    return compute_energy(ff, r), compute_energy_pullback
end

function optimize_structure!(ff::ForceField{T}; method=:OptimJL, kwargs...) where {T<:Real}
    if method == :OptimizationJL
        return optimize_structure_optimizationjl!(ff; kwargs...)
    elseif method == :OptimJL
        return optimize_structure_optimjl!(ff; kwargs...)
    end
end

function optimize_structure_optimizationjl!(
        ff::ForceField{T}; 
        algorithm=BFGS(), 
        kwargs...) where {T<:Real}
    r_flat = collect(Iterators.flatten(atoms_df(ff.system).r))

    optf = OptimizationFunction((r,_) -> compute_energy(ff, r), Optimization.AutoZygote())
    prob = OptimizationProblem(optf, r_flat)

    solve(prob, algorithm; kwargs...)
end

function optimize_structure_optimjl!(
        ff::ForceField{T}; 
        algorithm=BFGS(),
        kwargs...) where {T<:Real}
    r_flat = collect(Iterators.flatten(atoms_df(ff.system).r))

    optimize(
        r -> compute_energy(ff, r), 
        r -> Zygote.gradient(compute_energy, ff, r)[2], 
        r_flat, 
        method=algorithm, kwargs...)
end


function optimize_hydrogen_positions!(ff::ForceField{T}; method=:OptimJL, kwargs...) where {T<:Real}
    old_constrained = copy(ff.constrained_atoms)
    empty!(ff.constrained_atoms)
    append!(ff.constrained_atoms, findall(a -> a.element != Elements.H, atoms(ff.system)))
    solution = optimize_structure!(ff; method=method, kwargs...)
    empty!(ff.constrained_atoms)
    append!(ff.constrained_atoms, old_constrained)

    solution
end

function optimize_structure!(ff::Observable{ForceField{T}}; method=:OptimJL, notification_frequency=1, kwargs...) where {T<:Real}
    iterations = 0

    _update = () -> (
        iterations += 1;
        if iterations % notification_frequency == 0
            notify(ff)
        end;

        false
    )

    cb = if method == :OptimJL
        x -> _update()
    elseif method == :OptimizationJL
        (θ, l) -> _update()
    end

    optimize_structure!(ff[]; method=method, callback=cb, kwargs...)
end