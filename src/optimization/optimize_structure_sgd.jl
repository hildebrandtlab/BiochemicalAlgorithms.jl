

export
    optimize_structure_mini!

"""
    optimize_structure!(ff::ForceField)

Attempts to solve the energy optimization problem represented by the given force field object with a minibatching approach.

# Supported keyword arguments
This function passes all keyword arguments to
[Optimization.solve](https://docs.sciml.ai/Optimization/stable/API/solve/),
with the following default values:
 - `alg = ()`
"""
function optimize_structure_mini!(ff::ForceField; alg=OptimizationOptimisers.Adam(0.01), epochs::Int=10, batchsize::Int=10, callback=nothing, kwargs...)
    r0 = collect(Float64, Iterators.flatten(atoms(ff.system).r))

    ds = InteractionDataSet(ff)
    dataloader = MLUtils.DataLoader(ds, batchsize=batchsize, shuffle=true)
    batches = collect(dataloader)

    state = MiniBatchParams(ff, batches, 1)

    optf = Optimization.OptimizationFunction(
        (r, p=nothing) -> begin
            p = p !== nothing ? p : state
            _compute_energy_loss(r, p)
        end,
        grad = (g, r, p=nothing) -> begin
            p = p !== nothing ? p : state
            _compute_grad!(g, r, p)
        end
    )

    prob = Optimization.OptimizationProblem(optf, r0, state)

    iters_in_epoch = Ref(0)
    epoch_steps = Ref(max(1, length(state.batches)))
    done_epochs = Ref(0)

    combined_callback = (opt_state, l) -> begin
        stop = _epoch_minibatch_callback(
            opt_state, l, state, iters_in_epoch, epoch_steps, done_epochs, epochs, batchsize
        )
        if callback !== nothing
            e = compute_energy!(ff)
            callback(opt_state.iter, e)
        end
        opt_state.iter % 1000 == 0 && 
        return stop
    end

    sol = Optimization.solve(
        prob,
        alg;
        callback = combined_callback,
        maxiters = max(1, epochs * epoch_steps[]),
        kwargs...
    )

    return sol
end


struct InteractionBatch
    stretch_range::UnitRange{Int64}, 
    bend_range::UnitRange{Int64},
    torsion_range::UnitRange{Int64},
    improper_range::UnitRange{Int64}, 
    ljp_range::UnitRange{Int64}, 
    es_range::UnitRange{Int64}, 
    hb_range::UnitRange{Int64}
end


struct InteractionDataSet{T}
    ff::ForceField{T}
    number_atom_pairs::Int64
    data::Vector{InteractionBatch}

    function InteractionDataSet{T}(ff::ForceField{T}, batchsize) where T

        


        number_atom_pairs = 0

        data = Vector{InteractionBatch}()

        
        number_atom_pairs = length(data)
        new{T}(ff, number_atom_pairs, data)
    end
end

InteractionDataSet(ff::ForceField{T}) where T = InteractionDataSet{T}(ff)

Base.length(d::InteractionDataSet) = d.number_atom_pairs
Base.getindex(d::InteractionDataSet, i::Int) = d.data[i]
Base.getindex(d::InteractionDataSet, idxs::AbstractVector) = [getindex(d, i) for i in idxs]
