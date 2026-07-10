export
    InteractionBatch,
    InteractionDataSet,
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
function optimize_structure_mini!(cff::CompiledForceField{T,Acc}; alg=OptimizationOptimisers.Adam(0.01), epochs::Int=10, batchsize::Int=10, callback=nothing, kwargs...) where {T,Acc}
    r0 = collect(Float64, Iterators.flatten(atoms(cff.ff.system).r))

    ds = InteractionDataSet(cff, epochs, batchsize)

    state = MiniBatchParams(cff, ds.batches, 1)

    optf = Optimization.OptimizationFunction(
        (r, p=nothing) -> begin
            p = p !== nothing ? p : state
            _compute_energy_loss(r, p)
        end,
        grad=(g, r, p=nothing) -> begin
            p = p !== nothing ? p : state
            _compute_grad!(g, r, p)
        end
    )

    prob = Optimization.OptimizationProblem(optf, r0, state)

    iters_in_epoch = Ref(0)
    epoch_steps = Ref(max(1, length(ds.batches)))
    done_epochs = Ref(0)
    last_r = Ref(copy(r0))  # ← Track last parameters seen

    combined_callback = (opt_state, l) -> begin
        # Advance batch
        state.current_batch_idx = mod1(state.current_batch_idx + 1, length(state.batches))

        stop = _epoch_minibatch_callback(
            opt_state, l, state, iters_in_epoch, epoch_steps, done_epochs, epochs, batchsize
        )
        if callback !== nothing
            e = compute_energy!(cff)
            callback(opt_state.iter, e)
        end
        return stop
    end

    sol = Optimization.solve(
        prob,
        alg;
        callback=combined_callback,
        maxiters=max(1, epochs * epoch_steps[]),
        kwargs...
    )
    # write the final optimized state back into the canonical atom table
    finalize_optimization!(cff, sol.u)

    sol
end


struct InteractionBatch
    stretch_range::UnitRange{Int}
    bend_range::UnitRange{Int}
    torsion_range::UnitRange{Int}
    improper_range::UnitRange{Int}
    ljp_range::UnitRange{Int}
    es_range::UnitRange{Int}
    hb_range::UnitRange{Int}
end

function _partition_range(count::Integer, nbatches::Integer, batch::Integer)
    count <= 0 && return 1:0
    nbatches <= 1 && return 1:count

    base = div(count, nbatches)
    rem = count % nbatches
    start = 1
    for b in 1:nbatches
        size = base + (b <= rem ? 1 : 0)
        if b == batch
            return start:(start+size-1)
        end
        start += size
    end
    return 1:0
end

struct InteractionDataSet{T,Acc}
    cff::CompiledForceField{T,Acc}
    batchsize::Int
    epochs::Int
    batches::Vector{InteractionBatch}
    ninteractions::Int
    counts
    function InteractionDataSet{T,Acc}(cff::CompiledForceField{T,Acc}, batchsize::Integer=10, epochs::Integer=100) where {T,Acc}

        # generate a new permutation of the interactions
        shuffle_interactions!(cff)

        counts = (
            stretch=length(cff.stretch.i),
            bend=length(cff.bend.i),
            torsion=length(cff.proper.i),
            improper=length(cff.improper.i),
            ljp=length(cff.lj.i),
            es=length(cff.es.i),
            hb=length(cff.hb.i),
        )


        ninteractions = sum(counts)
        nbatches = max(1, cld(ninteractions, batchsize))
        batches = Vector{InteractionBatch}(undef, nbatches)


        for b in 1:nbatches
            batches[b] = InteractionBatch(
                _partition_range(counts.stretch, nbatches, b),
                _partition_range(counts.bend, nbatches, b),
                _partition_range(counts.torsion, nbatches, b),
                _partition_range(counts.improper, nbatches, b),
                _partition_range(counts.ljp, nbatches, b),
                _partition_range(counts.es, nbatches, b),
                _partition_range(counts.hb, nbatches, b),
            )
        end

        new{T,Acc}(cff, Int(batchsize), Int(epochs), batches, ninteractions, counts)
    end
end

InteractionDataSet(cff::CompiledForceField{T,Acc}, batchsize::Integer=10, epochs::Integer=100) where {T,Acc} =
    InteractionDataSet{T,Acc}(cff, batchsize, epochs)
InteractionDataSet(ff::ForceField{T}, batchsize::Integer=10, epochs::Integer=100) where {T} =
    InteractionDataSet(compile(ff), batchsize, epochs)

Base.length(d::InteractionDataSet) = d.ninteractions
Base.getindex(d::InteractionDataSet, i::Int) = d.batches[i]
Base.getindex(d::InteractionDataSet, idxs::AbstractVector) = [getindex(d, i) for i in idxs]


"""Mutable container for tracking batch index during minibatch optimization."""
mutable struct MiniBatchParams{T,Acc}
    cff::CompiledForceField{T,Acc}
    batches::Vector{InteractionBatch}
    current_batch_idx::Int
end

function _refresh_minibatches!(p::MiniBatchParams{T,Acc}, r0; batchsize::Int, shuffle::Bool=true) where {T,Acc}
    prepare_eval!(p.cff, r0)
    shuffle_interactions!(p.cff)
    ds = InteractionDataSet(p.cff, batchsize)
    p.batches = ds.batches
    p.current_batch_idx = 1
    return max(1, length(p.batches))
end

function _epoch_minibatch_callback(
    opt_state,
    l,
    p::MiniBatchParams{T,Acc},
    iters_in_epoch::Base.RefValue{Int},
    epoch_steps::Base.RefValue{Int},
    done_epochs::Base.RefValue{Int},
    epochs::Int,
    batchsize::Int
) where {T,Acc}
    iters_in_epoch[] += 1

    if iters_in_epoch[] >= epoch_steps[]
        done_epochs[] += 1
        done_epochs[] >= epochs && return true

        r0 = opt_state.u # pass the current state

        epoch_steps[] = _refresh_minibatches!(p, r0; batchsize=batchsize, shuffle=true)
        iters_in_epoch[] = 0
    end
    return false
end


function _compute_energy_loss(r::AbstractVector, p::MiniBatchParams{T,Acc}) where {T,Acc}
    prepare_eval!(p.cff, r)
    batch = p.batches[p.current_batch_idx]
    r = p.cff.r

    if p.cff.backend == BiochemicalAlgorithms.SerialBackend()
        stretchp = zero(Acc)
        bendp = zero(Acc)
        properp = zero(Acc)
        improperp = zero(Acc)
        vdwp = zero(Acc)
        hbp = zero(Acc)
        esp = zero(Acc)

        stretchp =  _sum_stretch(r, p.cff.stretch, first(batch.stretch_range), last(batch.stretch_range))
        bendp = _sum_bend(r, p.cff.bend, first(batch.bend_range), last(batch.bend_range))
        properp =_sum_torsion(r, p.cff.proper, first(batch.torsion_range), last(batch.torsion_range))
        improperp =  _sum_torsion(r, p.cff.improper, first(batch.improper_range), last(batch.improper_range))
        vdwp =_sum_lj(r, p.cff.lj, p.cff.vdw_sw, first(batch.ljp_range), last(batch.ljp_range))
        hbp =_sum_hb(r, p.cff.hb, p.cff.vdw_sw, first(batch.hb_range), last(batch.hb_range))
        esp =  _sum_es(r, p.cff.es, p.cff.es_sw, p.cff.distance_dependent_dielectric,
                      p.cff.es_prefactor, first(batch.es_range), last(batch.es_range))


        return sum(stretchp) + sum(bendp) + sum(properp) + sum(improperp) + sum(vdwp) + sum(hbp) + sum(esp)
    end

    nt = Threads.nthreads()
    stretchp = zeros(Acc, nt)
    bendp = zeros(Acc, nt)
    properp = zeros(Acc, nt)
    improperp = zeros(Acc, nt)
    vdwp = zeros(Acc, nt)
    hbp = zeros(Acc, nt)
    esp = zeros(Acc, nt)

    Threads.@threads :static for t in 1:nt
        sr = _chunk_range(batch.stretch_range, nt, t)
        
        stretchp[t] = _sum_stretch(r, p.cff.stretch, first(sr), last(sr))
    

        br = _chunk_range(batch.bend_range, nt, t)
        bendp[t] = _sum_bend(r, p.cff.bend, first(br), last(br))
   

        pr = _chunk_range(batch.torsion_range, nt, t)
        properp[t] = _sum_torsion(r, p.cff.proper, first(pr), last(pr))
       

        ir = _chunk_range(batch.improper_range, nt, t)
        
        improperp[t] = _sum_torsion(r, p.cff.improper, first(ir), last(ir))
       

        vr = _chunk_range(batch.ljp_range, nt, t)
        vdwp[t] = _sum_lj(r, p.cff.lj, p.cff.vdw_sw, first(vr), last(vr))
        
        hr = _chunk_range(batch.hb_range, nt, t)
        hbp[t] = _sum_hb(r, p.cff.hb, p.cff.vdw_sw, first(hr), last(hr))
        

        er = _chunk_range(batch.es_range, nt, t)
        esp[t] = _sum_es(r, p.cff.es, p.cff.es_sw, p.cff.distance_dependent_dielectric,
                             p.cff.es_prefactor, first(er), last(er))  
    end
   
   
    return  sum(stretchp) + sum(bendp) + sum(properp) + sum(improperp) + sum(vdwp) + sum(hbp) + sum(esp)
end

function _compute_grad!(grad::AbstractVector, r::AbstractVector, p::MiniBatchParams{T,Acc}) where {T,Acc}
    prepare_eval!(p.cff, r)
    batch = p.batches[p.current_batch_idx]

    F = p.cff.F
    fill!(F, zero(Vector3{Acc}))

    r = p.cff.r

   if p.cff.backend == BiochemicalAlgorithms.SerialBackend()
        _accum_stretch!(F, r, p.cff.stretch, batch.stretch_range.start, batch.stretch_range.stop)
        _accum_bend!(F, r, p.cff.bend, batch.bend_range.start, batch.bend_range.stop)
        _accum_torsion!(F, r, p.cff.proper, batch.torsion_range.start, batch.torsion_range.stop)
        _accum_torsion!(F, r, p.cff.improper, batch.improper_range.start, batch.improper_range.stop)
        _accum_lj!(F, r, p.cff.lj, p.cff.vdw_sw, batch.ljp_range.start, batch.ljp_range.stop)
        _accum_hb!(F, r, p.cff.hb, p.cff.vdw_sw, batch.hb_range.start, batch.hb_range.stop)
        _accum_es!(F, r, p.cff.es, p.cff.es_sw, p.cff.distance_dependent_dielectric,
            p.cff.es_prefactor_force,
            batch.es_range.start,
            batch.es_range.stop
        )
    else

        nt = Threads.nthreads()
            Threads.@threads :static for t in 1:nt
            
            Ft = p.cff.F_threads[t]
            fill!(Ft, zero(Vector3{Acc}))

            sr = _chunk_range(batch.stretch_range, nt, t)
            if !isempty(sr)
                _accum_stretch!(Ft, r, p.cff.stretch, first(sr), last(sr))
            end

            br = _chunk_range(batch.bend_range, nt, t)
            if !isempty(br)
                _accum_bend!(Ft, r, p.cff.bend, first(br), last(br))
            end

            pr = _chunk_range(batch.torsion_range, nt, t)
            if !isempty(pr)
                _accum_torsion!(Ft, r, p.cff.proper, first(pr), last(pr))
            end

            ir = _chunk_range(batch.improper_range, nt, t)
            if !isempty(ir)
                _accum_torsion!(Ft, r, p.cff.improper, first(ir), last(ir))
            end

            vr = _chunk_range(batch.ljp_range, nt, t)
            if !isempty(vr)
                _accum_lj!(Ft, r, p.cff.lj, p.cff.vdw_sw, first(vr), last(vr))
            end

            hr = _chunk_range(batch.hb_range, nt, t)
            if !isempty(hr)
                _accum_hb!(Ft, r, p.cff.hb, p.cff.vdw_sw, first(hr), last(hr))
            end

            er = _chunk_range(batch.es_range, nt, t)
            if !isempty(er)
                _accum_es!(Ft, r, p.cff.es, p.cff.es_sw, p.cff.distance_dependent_dielectric,
                        p.cff.es_prefactor_force, first(er), last(er))
            end

            @inbounds for t in 1:nt
                Ft = p.cff.F_threads[t]
                for k in 1:p.cff.natoms
                    F[k] += Ft[k]
                end
            end
        end
    end

    gradient_flat!(grad, p.cff)
    nothing
end

function _chunk_range(range::UnitRange{Int}, nt::Int, t::Int)
    isempty(range) && return first(range):(first(range) - 1)
    nt <= 1 && return range

    len = length(range)
    base = div(len, nt)
    rem = len % nt

    offset = (t - 1) * base + min(t - 1, rem)
    size = base + (t <= rem ? 1 : 0)

    start = first(range) + offset
    stop = start + size - 1
    start:stop
end
