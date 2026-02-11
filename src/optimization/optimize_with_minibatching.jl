export
    BatchParams,
    InteractionDataSet,
    getobs,
    _compute_energy_loss!,
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

struct BatchParams
    batch::Vector{IndexedInteraction}
    ff::ForceField
end 

function _compute_energy_loss!(r::Vector{Float64}, p::BatchParams)
    batch, ff = p.batch, p.ff        

    atoms(ff.system).r .= eachcol(reshape(r, 3, :))
    update!(ff)

    stretch_idx = [i[3] for i in batch if i[1] == Interaction.Stretch]
    bend_idx = [i[3] for i in batch if i[1] == Interaction.Bends]

    proper_torsions = [i[3] for i in batch if i[1] == Interaction.Torsion]
    improper_torsions = [i[3] for i in batch if i[1] == Interaction.ImproperTorsion]
    lj_interactions_idx = [i[3] for i in batch if i[1] == Interaction.LJP]
    hydrogen_bonds_idx = [i[3] for i in batch if i[1] == Interaction.HydrogenBond]
    electrostatic_interactions_idx = [i[3] for i in batch if i[1] == Interaction.Electrostatic]

    # TODO: call components to compute energy for batch
    map(idx -> compute_energy!(ff.components[1].stretches[idx]), stretch_idx)
    map(idx -> compute_energy!(ff.components[2].bends[idx]), bend_idx)
    map(compute_energy!, ff.components[2].bends[filter(i -> i[1] == Interaction.Bends, batch)])
    map(compute_energy!, ff.components[3].proper_torsions[filter(i -> i[1] == Interaction.Torsion, batch)])
    map(compute_energy!, ff.components[3].improper_torsions[filter(i -> i[1] == Interaction.ImproperTorsion, batch)])
    map(compute_energy!, [v for (i, v) in enumerate(ff.components[4].lj_interactions) if i in filter(j -> j[1] == Interaction.LJP, batch)])
    map(compute_energy!, [v for (i, v) in enumerate(ff.components[4].hydrogen_bonds) if i in filter(i -> i[1] == Interaction.HydrogenBonds, batch)])
    map(compute_energy!, [v for (i, v) in enumerate(ff.components[4].electrostatic_interactions) if i in filter(i -> i[1] == Interaction.Electrostatic, batch)])

    total_energy = 0

    for c in ff.components
        for (name, value) in c.energy
            ff.energy[name] = value
            total_energy += value
        end
    end

    if verbose
        @info "AMBER Energy:"

        max_length = maximum([length(k) for (k, _) in ff.energy])
        f_string = Printf.Format("%-$(max_length)s: %.9g kJ/mol")

        for (name, value) in ff.energy
            @info Printf.format(f_string, name, value)
        end

        @info repeat("-", max_length + 19)
        @info Printf.format(f_string, "total energy:", total_energy)
    end

    total_energy
end

function _compute_grad!(grad::Vector{Float64}, r::Vector{Float64}, p::BatchParams)
   batch, ff = p.batch, p.ff      
    update!(ff)

    # call components to compute forces
    map(compute_forces!, ff.components[1].stretches[filter(i -> i[1] == Interaction.Stretch, batch)])
    map(compute_forces!, ff.components[2].bends[filter(i -> i[1] == Interaction.Bends, batch)])


    map(compute_forces!, ff.components[3].proper_torsions[filter(i -> i[1] == Interaction.Torsion, batch)])
    map(compute_forces!, ff.components[3].improper_torsions[filter(i -> i[1] == Interaction.ImproperTorsion, batch)])


    map(compute_forces!, [v for (i, v) in enumerate(ff.components[4].lj_interactions) if i in filter(j -> j[1] == Interaction.LJP, batch)])
    map(compute_forces!, [v for (i, v) in enumerate(ff.components[4].hydrogen_bonds) if i in filter(i -> i[1] == Interaction.HydrogenBonds, batch)])
    map(compute_forces!, [v for (i, v) in enumerate(ff.components[4].electrostatic_interactions) if i in filter(i -> i[1] == Interaction.Electrostatic, batch)])


    # update forces
    F = atoms(ff.system).F
    F[ff.constrained_atoms] .= Ref(zeros(3))
    grad .= -collect(Float64, Iterators.flatten(F))
    nothing
end