export LennardJonesInteraction, ElectrostaticInteraction, NonBondedComponent

const ES_Prefactor = ustrip(Constants.N_A * Constants.e₀^2 / (4.0 * π * Constants.ε₀) |> u"angstrom*kJ/mol")

@auto_hash_equals struct LennardJonesInteraction{T<:Real, N, M}
    A::T
    B::T
    distance::T
    scaling_factor::T
    a1::Atom{T}
    a2::Atom{T}
    switching_function
end

@auto_hash_equals struct ElecrostaticInteraction{T<:Real}
    q1q2::T
    distance::T
    scaling_factor::T
    distance_dependent_dielectric::Bool
    a1::Atom{T}
    a2::Atom{T}
    switching_function
end

# the switching function is defined as follows:
#         (r_{off}^2 - R^2)^2 (r_{off}^2 + 2 R^2 - 3r_{on}^2)
# sw(R) = ---------------------------------------------------
#                    (r_{off}^2 - r_{on}^2)^3
#
# [Brooks et al., J. Comput. Chem., 4:191 (1983)]
#
# the derivative has the following form:
#
#                      (r_{off}^2 - R^2)(r_{on}^2 - R^2)
#   d/dR sw(R) = 12 R -----------------------------------
#                           (r_{off}^2 - r_{on}^2)^3
#
function _build_cubic_switching_function(cutoff::T, cuton::T) where {T<:Real}
    sq_cutoff = cutoff^2
    sq_cuton  = cuton^2

    inverse_distance_off_on_3 = (sq_cutoff - sq_cuton)^3

    if (inverse_distance_off_on_3 <= T(0.0))
        @warn "NonBondedComponent(): cuton value should be smaller than cutoff -- " *
                "switching function disabled."

        inverse_distance_off_on_3 = T(0.0)
    else
        inverse_distance_off_on_3 = T(1.0) / inverse_distance_off_on_3
    end

    # now, build the switching function as a closure
    sq_distance::T -> (
        below_off = ((sq_distance < sq_cutoff) ? one(T) : zero(T));
        below_on  = ((sq_distance < sq_cuton)  ? one(T) : zero(T));

        below_off * (below_on + (one(T) - below_on) * (sq_cutoff - sq_distance)^2
                  * (sq_cutoff + T(2.0) * sq_distance - T(3.0) * sq_cuton)
                  * inverse_distance_off_on_3)
    )
end

function _try_assign_vdw!(
        atom_1::Atom{T},
        atom_2::Atom{T},
        distance::T,
        scaling_factor::T,
        lj_combinations,
        lj_interactions,
        switching_function) where {T<:Real}

    type_atom_1 = atom_1.atom_type
    type_atom_2 = atom_2.atom_type

    params = get(lj_combinations, (I=type_atom_1, J=type_atom_2,), missing)

    if !ismissing(params)
        params = only(params)
        push!(
            lj_interactions,
                LennardJonesInteraction{T, 12, 6}(
                    params.A_ij,
                    params.B_ij,
                    T(distance),
                    scaling_factor,
                    atom_1,
                    atom_2,
                    switching_function
                )
        )
    else
        @warn "NonBondedComponent(): cannot find vdW parameters for "                   *
                "atom types $(type_atom_1)-$(type_atom_2) (atoms are: "                 *
                "$(get_full_name(atom_1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                "$(get_full_name(atom_2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"

        push!(ff.unassigned_atoms, atom_1)
        push!(ff.unassigned_atoms, atom_2)

        if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
            throw(TooManyErrors())
        end
    end
end

@auto_hash_equals struct NonBondedComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    lj_interactions::AbstractVector{LennardJonesInteraction{T, 12, 6}}
    hydrogen_bonds::AbstractVector{LennardJonesInteraction{T, 12, 10}}
    electrostatic_interactions::AbstractVector{ElecrostaticInteraction{T}}
    vdw_switching_function
    es_switching_function
    num_1_4_interactions::T
    energy::Dict{String, T}

    function NonBondedComponent{T}(ff::ForceField{T}) where {T<:Real}
        # extract the parameter section for quadratic angle bends
        lj_section = extract_section(ff.parameters, "LennardJones")
        lj_df = lj_section.data

        unit_R = get(lj_section.properties, "unit_R",       "angstrom")
        unit_ϵ = get(lj_section.properties, "unit_epsilon", "kcal/mol")

        # ball used to write Angstrom with a capital letter; this clashes with the convention in Unitful.jl
        if unit_R == "Angstrom"
            unit_R = "angstrom"
        end

        R_factor = ustrip((1uparse(unit_R)) |> u"angstrom")
        ϵ_factor = ustrip((1uparse(unit_ϵ)) |> u"kJ/mol")

        lj_mininmal_df = lj_df[:, [:I, :R, :epsilon]]
        lj_mininmal_df.R *= R_factor
        lj_mininmal_df.epsilon *= ϵ_factor

        lj_combinations = crossjoin(lj_mininmal_df, lj_mininmal_df; makeunique=true)
        rename!(lj_combinations, 
            [:R => :R_i, :epsilon => :epsilon_i, :I_1 => :J, :R_1 => :R_j, :epsilon_1 => :epsilon_j])

        # convert the data from RE - format to A/B - format
        radius_averaging = get(lj_section.properties, :radius_averaging, "arithmetic")

        if radius_averaging ∉ ["arithmetic", "geometric"]
            @warn "NonBondedComponent(): unknown radius averaging method $(radius_averaging)."
            @warn "Using arithmetic averaging instead."

            radius_averaging = "arithmetic"
        end

        r_6 = ((radius_averaging == "arithmetic")
                ? (lj_combinations.R_i .+ lj_combinations.R_j)
             : T(2.0) * sqrt.(lj_combinations.R_i .* lj_combinations.Rj)
            ).^6

        ϵ = sqrt.(lj_combinations.epsilon_i .* lj_combinations.epsilon_j)

        lj_combinations.A_ij = ϵ .* r_6.^2
        lj_combinations.B_ij = T(2.0) * ϵ .* r_6
        
        lj_combinations = groupby(lj_combinations, ["I", "J"])

        # we will also need the hydrogen bond section to determine if parameters are missing for a pair of atoms
        hydrogen_bonds_section = extract_section(ff.parameters, "HydrogenBonds")
        hydrogen_bonds_df = hydrogen_bonds_section.data

        unit_hbond_A = get(hydrogen_bonds_section.properties, "unit_A", "kcal/mol*angstrom^12")
        unit_hbond_B = get(hydrogen_bonds_section.properties, "unit_B", "kcal/mol*angstrom^10")

        # here, BALL used A for angstrom; this clashes with the convention in Unitful.jl
        unit_hbond_A = replace(unit_hbond_A, "A" => "angstrom")
        unit_hbond_B = replace(unit_hbond_B, "A" => "angstrom")

        hbond_A_factor = ustrip((1uparse(unit_hbond_A)) |> u"kJ/mol*angstrom^12")
        hbond_B_factor = ustrip((1uparse(unit_hbond_B)) |> u"kJ/mol*angstrom^10")

        # hbond terms are already stored in A/B - format, so conversion is simple
        hydrogen_bonds_df.A .*= T(hbond_A_factor)
        hydrogen_bonds_df.B .*= T(hbond_B_factor)

        # group the hydrogen bond parameters by type_i, type_j combinations
        hydrogen_bond_combinations = groupby(hydrogen_bonds_df, ["I", "J"])

        # the cutoffs for the nonbonded pair list and the switching function
        nonbonded_cutoff = ff.options[:nonbonded_cutoff]
        vdw_cutoff       = ff.options[:vdw_cutoff]
        vdw_cuton        = ff.options[:vdw_cuton]
        es_cutoff        = ff.options[:electrostatic_cutoff]
        es_cuton         = ff.options[:electrostatic_cuton]  

        # build the switching function as a closure
        vdw_switching_function = _build_cubic_switching_function(vdw_cutoff, vdw_cuton)
        es_switching_function  = _build_cubic_switching_function(es_cutoff, es_cuton)
        
        scaling_vdw_1_4 = ff.options[:scaling_vdw_1_4]
        if (scaling_vdw_1_4 == T(0.0))
            @warn "NonBondedComponent(): illegal - 1-4 vdW scaling factor: must be non-zero!"
            @warn "Resetting to 1.0."

            scaling_vdw_1_4 = T(1.0)
        else
            scaling_vdw_1_4 = T(1.0) / scaling_vdw_1_4
        end

        scaling_es_1_4 = ff.options[:scaling_electrostatic_1_4]
        if (scaling_es_1_4 == T(0.0))
            @warn "NonBondedComponent(): illegal - 1-4 electrostatic scaling factor: must be non-zero!"
            @warn "Resetting to 1.0."

            scaling_es_1_4 = T(1.0)
        else
            scaling_es_1_4 = T(1.0) / scaling_es_1_4
        end

        # do we need to construct a periodic box?
        periodic_box = [
            ff.options[:periodic_box_width], 
            ff.options[:periodic_box_height],
            ff.options[:periodic_box_depth]
        ]

        neighbors = ((ff.options[:periodic_boundary_conditions]) 
            ? neighborlist(atoms_df(ff.system).r, unitcell=periodic_box, nonbonded_cutoff)
            : neighborlist(atoms_df(ff.system).r, nonbonded_cutoff)
        )

        distance_dependent_dielectric = ff.options[:distance_dependent_dielectric]


        lj_interactions = Vector{LennardJonesInteraction{T, 12, 6}}()
        hydrogen_bonds  = Vector{LennardJonesInteraction{T, 12, 10}}()
        electrostatic_interactions = Vector{ElecrostaticInteraction{T}}()

		# count the number of 1-4 interactions (torsions)
		num_1_4_interactions = 0

        for lj_candidate in neighbors
            atom_1 = atoms(ff.system)[lj_candidate[1]]
            atom_2 = atoms(ff.system)[lj_candidate[2]]

            # exclude 1-2 and 1-3 interactions
            if is_bound_to(atom_1, atom_2) || is_geminal(atom_1, atom_2)
                continue
            end

            vicinal_pair = is_vicinal(atom_1, atom_2)

            q1q2 = atom_1.charge * atom_2.charge

            if q1q2 ≠ zero(T)
                push!(
                    electrostatic_interactions,
                    ElecrostaticInteraction{T}(
                        atom_1.charge * atom_2.charge,
                        T(lj_candidate[3]),
                        vicinal_pair ? scaling_es_1_4 : T(1.0),
                        distance_dependent_dielectric,
                        atom_1,
                        atom_2,
                        es_switching_function
                    )
                )
            end

            # first, figure out if the atoms are part of a torsion
            if !vicinal_pair              
                # no. now, figure out if these atoms are part of a hydrogen bond

                # parameters for hydrogen bonds are used, if they exist
                # and the two atoms are not vicinal (1-4).

                h_params = coalesce(
                    get(hydrogen_bond_combinations, (I=atom_1.atom_type, J=atom_2.atom_type,), missing),
                    get(hydrogen_bond_combinations, (I=atom_2.atom_type, J=atom_1.atom_type,), missing)
                )

                if !ismissing(h_params)
                    params = only(params)
                    push!(
                        hydrogen_bonds,
                        LennardJonesInteraction{T, 12, 10}(
                            params.A,
                            params.B,
                            T(lj_candidate[3]),
                            T(1.0),
                            atom_1,
                            atom_2,
                            vdw_switching_function
                        )
                    )
                else
                    _try_assign_vdw!(
                        atom_1,
                        atom_2,
                        T(lj_candidate[3]),
                        T(1.0),
                        lj_combinations,
                        lj_interactions,
                        vdw_switching_function
                    )
                end
            else
                # this is a torsion
                num_1_4_interactions += 1
                _try_assign_vdw!(
                    atom_1,
                    atom_2,
                    T(lj_candidate[3]),
                    scaling_vdw_1_4,
                    lj_combinations,
                    lj_interactions,
                    vdw_switching_function
                )
            end
        end

        new(
            "NonBonded", 
            ff,
            lj_interactions,
            hydrogen_bonds,
            electrostatic_interactions,
            vdw_switching_function,
            es_switching_function,
            num_1_4_interactions, 
            Dict{String, T}()
        )
    end
end

@inline function compute_energy(lji::LennardJonesInteraction{T, 12, 6}) where {T<:Real}
    inv_dist_6 = lji.distance^-6

    inv_dist_6 * (inv_dist_6 * lji.A - lji.B) * lji.scaling_factor * lji.switching_function(lji.distance^2)
end

@inline function compute_energy(hb::LennardJonesInteraction{T, 12, 10}) where {T<:Real}
    hb.distance^-12 * hb.A - hb.distance^-10 * hb.B * hb.scaling_factor * hb.switching_function(hb.distance^2)
end

@inline function compute_energy(esi::ElecrostaticInteraction{T}) where {T<:Real}
    energy = esi.distance_dependent_dielectric ? esi.q1q2 * esi.distance^-2 : esi.q1q2 / esi.distance

    energy * esi.scaling_factor * esi.switching_function(esi.distance^2) * T(ES_Prefactor)
end

function compute_energy(nbc::NonBondedComponent{T}) where {T<:Real}
    # iterate over all interactions in the system
    vdw_energy   = mapreduce(compute_energy, +, nbc.lj_interactions;            init=zero(T))
    hbond_energy = mapreduce(compute_energy, +, nbc.hydrogen_bonds;             init=zero(T))
    es_energy    = mapreduce(compute_energy, +, nbc.electrostatic_interactions; init=zero(T))

    nbc.energy["Van der Waals"]  = vdw_energy
    nbc.energy["Hydrogen Bonds"] = hbond_energy
    nbc.energy["Electrostatic"]  = es_energy
    
    vdw_energy + hbond_energy + es_energy
end

