export 
    InteractionDataSet


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
const IndexedInteraction = Tuple{InteractionType, AtomPair, Int} # Int is for the index of the corresponding interaction in the force field component

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
            push!(data, (Interaction.LJP, ap,  i)) 
        end
        #check for hydrogen bonds
        for (i, ljp) in enumerate(ff.components[4].hydrogen_bonds)
            ap = ljp.a1.idx < ljp.a2.idx ? (ljp.a1.idx, ljp.a2.idx) : (ljp.a2.idx, ljp.a1.idx)
            push!(data, (Interaction.HydrogenBond, ap, i))   
        end
        #check for electrostatic interactions
        for (i, ljp) in enumerate(ff.components[4].electrostatic_interactions)
            ap = ljp.a1.idx < ljp.a2.idx ? (ljp.a1.idx, ljp.a2.idx) : (ljp.a2.idx, ljp.a1.idx)
            push!(data, (Interaction.Electrostatic, ap,i))
        end     

        number_atom_pairs = length(data)
        new{T}(ff, number_atom_pairs, data)
    end
end

InteractionDataSet(ff::ForceField{T}) where T = InteractionDataSet{T}(ff)

Base.length(d::InteractionDataSet) = d.number_atom_pairs
Base.getindex(d::InteractionDataSet, i::Int) = d.data[i]    

function compute_forces!(ff::ForceField{T}; idset::InteractionDataSet{T}) where T
    # first, zero out the current forces
    atoms(ff.system).F .= Ref(zero(Vector3{T}))

    map(compute_forces!, ff.components)

    nothing
end

function compute_energy!(ff::ForceField{T}; verbose::Bool=false)::T where T
    total_energy = mapreduce(compute_energy!, +, ff.components; init=zero(T))

    for c in ff.components
        for (name, value) in c.energy
            ff.energy[name] = value
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