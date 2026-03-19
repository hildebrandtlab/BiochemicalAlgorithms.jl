export
    MiniBatchParams,
    InteractionDataSet,
    getobs,
    _compute_energy_loss,
    _compute_grad!


@enumx Interaction begin
    Stretch = 1
    Bends = 2
    Torsion = 3
    ImproperTorsion = 4
    LJP = 5
    HydrogenBond = 6
    Electrostatic = 7
    Unknown = 100
end

const InteractionType = Interaction.T
const AtomPair = Tuple{Int,Int}         # (atom index 1, atom index 2)
const IndexedInteraction = Tuple{InteractionType,AtomPair,Int} # Int is for the index of the corresponding interaction in the force field component

struct InteractionDataSet{T}
    ff::ForceField{T}
    number_atom_pairs::Int64
    data::Vector{IndexedInteraction}

    function InteractionDataSet{T}(ff::ForceField{T}) where T

        number_atom_pairs = 0

        data = Vector{IndexedInteraction}()

        # check stretches
        for (i, stretch) in enumerate(ff.components[1].stretches)
            ap = stretch.a1.idx < stretch.a2.idx ? (stretch.a1.idx, stretch.a2.idx) : (stretch.a2.idx, stretch.a1.idx)
            push!(data, (Interaction.Stretch, ap, i))
        end

        # check bends
        for (i, bend) in enumerate(ff.components[2].bends)
            ap = bend.a1.idx < bend.a2.idx ? (bend.a1.idx, bend.a2.idx) : (bend.a2.idx, bend.a1.idx)
            push!(data, (Interaction.Bends, ap, i))
        end
        # check torsions
        for (i, torsion) in enumerate(ff.components[3].proper_torsions)
            ap = torsion.a1.idx < torsion.a2.idx ? (torsion.a1.idx, torsion.a2.idx) : (torsion.a2.idx, torsion.a1.idx)
            push!(data, (Interaction.Torsion, ap, i))
        end
        # check improper torsions
        for (i, torsion) in enumerate(ff.components[3].improper_torsions)
            ap = torsion.a1.idx < torsion.a2.idx ? (torsion.a1.idx, torsion.a2.idx) : (torsion.a2.idx, torsion.a1.idx)
            push!(data, (Interaction.ImproperTorsion, ap, i))
        end
        #check for ljp
        for (i, ljp) in enumerate(ff.components[4].lj_interactions)
            ap = ljp.a1.idx < ljp.a2.idx ? (ljp.a1.idx, ljp.a2.idx) : (ljp.a2.idx, ljp.a1.idx)
            push!(data, (Interaction.LJP, ap, i))
        end
        #check for hydrogen bonds
        for (i, ljp) in enumerate(ff.components[4].hydrogen_bonds)
            ap = ljp.a1.idx < ljp.a2.idx ? (ljp.a1.idx, ljp.a2.idx) : (ljp.a2.idx, ljp.a1.idx)
            push!(data, (Interaction.HydrogenBond, ap, i))
        end
        #check for electrostatic interactions
        for (i, ljp) in enumerate(ff.components[4].electrostatic_interactions)
            ap = ljp.a1.idx < ljp.a2.idx ? (ljp.a1.idx, ljp.a2.idx) : (ljp.a2.idx, ljp.a1.idx)
            push!(data, (Interaction.Electrostatic, ap, i))
        end

        number_atom_pairs = length(data)
        new{T}(ff, number_atom_pairs, data)
    end
end

InteractionDataSet(ff::ForceField{T}) where T = InteractionDataSet{T}(ff)

Base.length(d::InteractionDataSet) = d.number_atom_pairs
Base.getindex(d::InteractionDataSet, i::Int) = d.data[i]
Base.getindex(d::InteractionDataSet, idxs::AbstractVector) = [getindex(d, i) for i in idxs]

function idx_by_type(batch::AbstractVector)
    return Dict(
        itype => [i[3] for i in batch if i[1] == itype]
        for itype in unique(i[1] for i in batch)
    )
end

"""Mutable container for tracking batch index during minibatch optimization."""
mutable struct MiniBatchParams
    ff::ForceField
    batches::Vector
    current_batch_idx::Int
end 

function _compute_energy_loss(r::Vector{T}, p::MiniBatchParams) where T
    # Cycle to next batch
    p.current_batch_idx = mod1(p.current_batch_idx + 1, length(p.batches))

    batch = p.batches[p.current_batch_idx]
    ff = p.ff

    atoms(ff.system).r .= eachcol(reshape(r, 3, :))

    # Build index mapping once
    idx_map = idx_by_type(batch)
    
    # Compute energy ONLY for batch interactions and update component energies
    total_energy = zero(T)
    
    # Stretches
    e_stretch = mapreduce(i -> compute_energy(ff.components[1].stretches[i]), +, get(idx_map, Interaction.Stretch, []); init=zero(T))
    total_energy += e_stretch
    
    # Bends
    e_bend = mapreduce(i -> compute_energy(ff.components[2].bends[i]), +, get(idx_map, Interaction.Bends, []); init=zero(T))
    total_energy += e_bend
    
    # Proper torsions
    e_proper = mapreduce(i -> compute_energy(ff.components[3].proper_torsions[i]), +, get(idx_map, Interaction.Torsion, []); init=zero(T))
    total_energy += e_proper
    
    # Improper torsions
    e_improper = mapreduce(i -> compute_energy(ff.components[3].improper_torsions[i]), +, get(idx_map, Interaction.ImproperTorsion, []); init=zero(T))
    total_energy += e_improper
    
    # Non-bonded interactions - enumerate and filter
    e_lj = mapreduce(
        v -> compute_energy(v), +,
        [v for (i, v) in enumerate(ff.components[4].lj_interactions) if i in get(idx_map, Interaction.LJP, [])];
        init=zero(T)
    )
    
    e_hbond = mapreduce(
        v -> compute_energy(v), +,
        [v for (i, v) in enumerate(ff.components[4].hydrogen_bonds) if i in get(idx_map, Interaction.HydrogenBond, [])];
        init=zero(T)
    )
    
    e_es = mapreduce(
        v -> compute_energy(v), +,
        [v for (i, v) in enumerate(ff.components[4].electrostatic_interactions) if i in get(idx_map, Interaction.Electrostatic, [])];
        init=zero(T)
    )
    total_energy += e_lj + e_hbond + e_es

    total_energy
end

function _compute_grad!(grad::Vector{T}, r::Vector{T}, p::MiniBatchParams) where T
    batch = p.batches[p.current_batch_idx]
    ff = p.ff
    
    atoms(ff.system).r .= eachcol(reshape(r, 3, :))
 #   update!(ff)

    # Build index mapping once
    idx_map = idx_by_type(batch)
    
    # Zero out forces before computing batch-selective forces
    F = atoms(ff.system).F
    F .= Ref(zeros(3))
    
    # Compute forces for batch interactions only
    f_stretch = map(compute_forces!, ff.components[1].stretches[get(idx_map, Interaction.Stretch, [])])
    f_bend = map(compute_forces!, ff.components[2].bends[get(idx_map, Interaction.Bends, [])])
    f_proper = map(compute_forces!, ff.components[3].proper_torsions[get(idx_map, Interaction.Torsion, [])])
    f_improper = map(compute_forces!, ff.components[3].improper_torsions[get(idx_map, Interaction.ImproperTorsion, [])])
    
    # Non-bonded interactions - iterate directly over filtered interactions
    for v in [v for (i, v) in enumerate(ff.components[4].lj_interactions) if i in get(idx_map, Interaction.LJP, [])]
        compute_forces!(v)
    end
    for v in [v for (i, v) in enumerate(ff.components[4].hydrogen_bonds) if i in get(idx_map, Interaction.HydrogenBond, [])]
        compute_forces!(v)
    end
    for v in [v for (i, v) in enumerate(ff.components[4].electrostatic_interactions) if i in get(idx_map, Interaction.Electrostatic, [])]
        compute_forces!(v)
    end

    # Extract gradients from batch forces
 #   F[ff.constrained_atoms] .= Ref(zeros(3))
    grad .= -collect(Float64, Iterators.flatten(F))
    nothing
end