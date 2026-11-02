export
    EEF1Interaction,
    CharmmNonBondedComponent

# EEF1 implicit solvation (Lazaridis & Karplus, Proteins 35:133 (1999)):
#
#   ΔG_solv = Σ_i ΔG_i^ref - Σ_i Σ_{j≠i} f_i(r_ij) V_j
#
# with the Gaussian solvation free energy density
#
#                2 ΔG_i^free
#   f_i(r) = -------------------- exp(-((r - R_i^min) / σ_i)^2)
#             4 π^(3/2) σ_i r^2
#
# Each EEF1Interaction stores the pair term f_i(r_ij) V_j + f_j(r_ij) V_i.
@auto_hash_equals struct EEF1Interaction{T<:Real}
    V1::T
    V2::T
    ΔG_free1::T
    ΔG_free2::T
    σ1::T
    σ2::T
    R_min1::T
    R_min2::T
    distance::T
    a1::Atom{T}
    a2::Atom{T}
    switching_function::CubicSwitchingFunction{T}
end

@auto_hash_equals mutable struct CharmmNonBondedComponentCache{T}
    bond_cache::Set{Pair{Int64, Int64}}
    es_switching_function::CubicSwitchingFunction{T}
    eef1_parameters::Dict{String, @NamedTuple{V::T, ΔG_ref::T, ΔG_free::T, σ::T, R_min::T}}
    geminal_cache::Set{Pair{Int64, Int64}}
    lj_combinations::Dict{@NamedTuple{I::String, J::String}, @NamedTuple{A_ij::T, B_ij::T}}
    lj_14_combinations::Dict{@NamedTuple{I::String, J::String}, @NamedTuple{A_ij::T, B_ij::T}}
    nonbonded_cutoff::T
    periodic_box::Vector{T}
    scaling_es_1_4::T
    scaling_vdw_1_4::T
    solvation_switching_function::CubicSwitchingFunction{T}
    use_eef1::Bool
    vdw_switching_function::CubicSwitchingFunction{T}
    vicinal_cache::Set{Pair{Int64, Int64}}

    function CharmmNonBondedComponentCache{T}() where T
        new()
    end
end

@auto_hash_equals mutable struct CharmmNonBondedComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    cache::CharmmNonBondedComponentCache{T}
    energy::Dict{String, T}
    unassigned_lj_interactions::Vector{Tuple{Atom{T}, Atom{T}}}

    solvation_reference_energy::T

    lj_interactions::Deque{LennardJonesInteraction{T, 12, 6}}
    electrostatic_interactions::Deque{ElectrostaticInteraction{T}}
    solvation_interactions::Deque{EEF1Interaction{T}}

    function CharmmNonBondedComponent{T}(ff::ForceField{T}) where {T<:Real}
        new("CharmmNonBonded", ff, CharmmNonBondedComponentCache{T}(), Dict{String, T}(), Tuple{Atom{T},Atom{T}}[], zero(T))
    end
end

# CHARMM omits the vdW interactions between the diagonal ring atoms of PHE and TYR
const _CHARMM_RING_DIAGONALS = Set([
    "CG"  => "CZ",  "CZ"  => "CG",
    "CD1" => "CE2", "CE2" => "CD1",
    "CD2" => "CE1", "CE1" => "CD2"
])

function _is_ring_diagonal(atom_1::Atom, atom_2::Atom)
    (atom_1.name => atom_2.name) ∈ _CHARMM_RING_DIAGONALS || return false

    frag = parent_fragment(atom_1)

    !isnothing(frag) &&
        atom_2.fragment_idx == frag.idx &&
        frag.name ∈ ("PHE", "TYR")
end

function _setup_eef1!(nbc::CharmmNonBondedComponent{T}) where T
    ff = nbc.ff

    nbc.cache.use_eef1 = ff.options[:use_eef1]::Bool && haskey(ff.parameters.sections, "EEF1Solvation")
    nbc.cache.eef1_parameters = Dict{String, @NamedTuple{V::T, ΔG_ref::T, ΔG_free::T, σ::T, R_min::T}}()
    nbc.solvation_reference_energy = zero(T)

    if !nbc.cache.use_eef1
        return
    end

    eef1_section = extract_section(ff.parameters, "EEF1Solvation")

    # ball used to write Angstrom with a capital letter; this clashes with the convention in Unitful.jl
    _unit_factor(name, default, target) = T(ustrip(
        (1uparse(replace(get(eef1_section.properties, name, default), "Angstrom" => "angstrom"))) |> target
    ))

    V_factor      = _unit_factor("unit_V",       "angstrom^3", u"angstrom^3")
    ΔG_ref_factor = _unit_factor("unit_dG_ref",  "kJ/mol",     u"kJ/mol")
    ΔG_free_factor= _unit_factor("unit_dG_free", "kJ/mol",     u"kJ/mol")
    σ_factor      = _unit_factor("unit_sig_w",   "angstrom",   u"angstrom")
    R_min_factor  = _unit_factor("unit_R_min",   "angstrom",   u"angstrom")

    for r in eachrow(eef1_section.data)
        nbc.cache.eef1_parameters[String(r.type)] = (
            V       = V_factor       * T(r.V),
            ΔG_ref  = ΔG_ref_factor  * T(r.dG_ref),
            ΔG_free = ΔG_free_factor * T(r.dG_free),
            σ       = σ_factor       * T(r.sig_w),
            R_min   = R_min_factor   * T(r.R_min)
        )
    end

    # the reference solvation energy only depends on the atom types
    for atom in atoms(ff.system)
        p = get(nbc.cache.eef1_parameters, atom.atom_type, missing)

        if ismissing(p)
            push!(ff.warnings,
                "CharmmNonBondedComponent(): no solvation parameters for atom type " *
                "$(atom.atom_type) ($(get_full_name(atom)))"
            )
        else
            nbc.solvation_reference_energy += p.ΔG_ref
        end
    end
end

function setup!(nbc::CharmmNonBondedComponent{T}) where T
    ff = nbc.ff

    atom_types = String.(extract_section(ff.parameters, "AtomTypes").data.type)

    nbc.cache.lj_combinations = _lj_combinations(
        extract_section(ff.parameters, "LennardJones"), T; types = atom_types
    )

    nbc.cache.lj_14_combinations = haskey(ff.parameters.sections, "LennardJones14") ?
        _lj_combinations(extract_section(ff.parameters, "LennardJones14"), T; types = atom_types) :
        empty(nbc.cache.lj_combinations)

    nbc.cache.vdw_switching_function = CubicSwitchingFunction{T}(
        ff.options[:vdw_cutoff]::T, ff.options[:vdw_cuton]::T
    )
    nbc.cache.es_switching_function = CubicSwitchingFunction{T}(
        ff.options[:electrostatic_cutoff]::T, ff.options[:electrostatic_cuton]::T
    )
    nbc.cache.solvation_switching_function = CubicSwitchingFunction{T}(
        ff.options[:solvation_cutoff]::T, ff.options[:solvation_cuton]::T
    )

    # in contrast to AMBER, CHARMM multiplies the 1-4 interactions with the scaling factors
    nbc.cache.scaling_es_1_4  = ff.options[:scaling_electrostatic_1_4]::T
    nbc.cache.scaling_vdw_1_4 = ff.options[:scaling_vdw_1_4]::T

    _setup_eef1!(nbc)

    nbc.cache.nonbonded_cutoff = ff.options[:nonbonded_cutoff]::T
    nbc.cache.periodic_box = [
        ff.options[:periodic_box_width]::T,
        ff.options[:periodic_box_height]::T,
        ff.options[:periodic_box_depth]::T
    ]

    bond_cache, geminal_cache, vicinal_cache = _bonded_neighbourhood_caches(ff.system)

    nbc.cache.bond_cache    = bond_cache
    nbc.cache.geminal_cache = geminal_cache
    nbc.cache.vicinal_cache = vicinal_cache
end

function update!(nbc::CharmmNonBondedComponent{T}) where T
    ff = nbc.ff
    cache = nbc.cache

    check_bond(a1::Int, a2::Int)    = (a1 => a2) ∈ cache.bond_cache
    check_geminal(a1::Int, a2::Int) = (a1 => a2) ∈ cache.geminal_cache
    check_vicinal(a1::Int, a2::Int) = (a1 => a2) ∈ cache.vicinal_cache

    atom_cache::AtomTable{T} = atoms(ff.system)
    neighbors::Vector{Tuple{Int, Int, T}} = ff.options[:periodic_boundary_conditions]::Bool ?
        neighborlist(positions=atom_cache.r, cutoff=cache.nonbonded_cutoff, unitcell=cache.periodic_box) :
        neighborlist(positions=atom_cache.r, cutoff=cache.nonbonded_cutoff)

    distance_dependent_dielectric = ff.options[:distance_dependent_dielectric]::Bool

    # ElectrostaticInteraction uses AMBER's distance dependent dielectric (epsilon = 4r),
    # while CHARMM uses epsilon = r; we compensate for this by scaling the charge product
    q1q2_factor = distance_dependent_dielectric ? T(4) : one(T)

    hint = Int(round(1.2 * natoms(ff.system)))
    lj_interactions = Deque{LennardJonesInteraction{T, 12, 6}}(hint)
    electrostatic_interactions = Deque{ElectrostaticInteraction{T}}(hint)
    solvation_interactions = Deque{EEF1Interaction{T}}(hint)
    empty!(nbc.unassigned_lj_interactions)

    for (i, j, d) in neighbors
        atom_1 = atom_cache[i]
        atom_2 = atom_cache[j]

        # exclude 1-2 and 1-3 interactions
        if check_bond(atom_1.idx, atom_2.idx) || check_geminal(atom_1.idx, atom_2.idx)
            continue
        end

        vicinal_pair = check_vicinal(atom_1.idx, atom_2.idx)
        distance = T(d)

        q1q2 = atom_1.charge * atom_2.charge

        if q1q2 ≠ zero(T)
            push!(electrostatic_interactions,
                ElectrostaticInteraction{T}(
                    q1q2_factor * q1q2,
                    distance,
                    vicinal_pair ? cache.scaling_es_1_4 : one(T),
                    distance_dependent_dielectric,
                    atom_1,
                    atom_2,
                    cache.es_switching_function
                )
            )
        end

        if vicinal_pair
            # 1-4 interactions use their own parameters if available
            if !_is_ring_diagonal(atom_1, atom_2)
                key = (I=atom_1.atom_type, J=atom_2.atom_type,)
                lj_combinations = haskey(cache.lj_14_combinations, key) ?
                    cache.lj_14_combinations :
                    cache.lj_combinations

                _try_assign_vdw!(
                    nbc,
                    atom_1,
                    atom_1.atom_type,
                    atom_2,
                    atom_2.atom_type,
                    distance,
                    cache.scaling_vdw_1_4,
                    lj_combinations,
                    lj_interactions,
                    cache.vdw_switching_function
                )
            end
        else
            _try_assign_vdw!(
                nbc,
                atom_1,
                atom_1.atom_type,
                atom_2,
                atom_2.atom_type,
                distance,
                one(T),
                cache.lj_combinations,
                lj_interactions,
                cache.vdw_switching_function
            )
        end

        # hydrogen atoms are not considered in EEF1
        if cache.use_eef1 && atom_1.element != Elements.H && atom_2.element != Elements.H
            p1 = get(cache.eef1_parameters, atom_1.atom_type, missing)
            p2 = get(cache.eef1_parameters, atom_2.atom_type, missing)

            if !ismissing(p1) && !ismissing(p2)
                push!(solvation_interactions,
                    EEF1Interaction{T}(
                        p1.V, p2.V,
                        p1.ΔG_free, p2.ΔG_free,
                        p1.σ, p2.σ,
                        p1.R_min, p2.R_min,
                        distance,
                        atom_1,
                        atom_2,
                        cache.solvation_switching_function
                    )
                )
            end
        end
    end

    nbc.lj_interactions            = lj_interactions
    nbc.electrostatic_interactions = electrostatic_interactions
    nbc.solvation_interactions     = solvation_interactions

    nothing
end

# returns the pair contribution to the solvation energy at the given distance
@inline function _eef1_energy(si::EEF1Interaction{T}, distance::T)::T where T
    inv_factor = one(T) / (T(π)^T(1.5) * distance^2)

    -T(0.5) * inv_factor * (
        si.V2 * si.ΔG_free1 * exp(-((distance - si.R_min1) / si.σ1)^2) / si.σ1 +
        si.V1 * si.ΔG_free2 * exp(-((distance - si.R_min2) / si.σ2)^2) / si.σ2
    )
end

@inline function compute_energy(si::EEF1Interaction{T})::T where T
    _eef1_energy(si, si.distance) * switching_function(si.switching_function, si.distance^2)
end

function compute_energy!(nbc::CharmmNonBondedComponent{T})::T where T
    vdw_energy = mapreduce(compute_energy, +, nbc.lj_interactions;            init=zero(T))
    es_energy  = mapreduce(compute_energy, +, nbc.electrostatic_interactions; init=zero(T))

    solvation_energy = nbc.cache.use_eef1 ?
        nbc.solvation_reference_energy + mapreduce(compute_energy, +, nbc.solvation_interactions; init=zero(T)) :
        zero(T)

    nbc.energy["Van der Waals"] = vdw_energy
    nbc.energy["Electrostatic"] = es_energy
    nbc.energy["Solvation"]     = solvation_energy

    vdw_energy + es_energy + solvation_energy
end

function compute_forces!(si::EEF1Interaction{T}) where T
    direction = si.a1.r .- si.a2.r

    sq_distance = squared_norm(direction)

    if (sq_distance > zero(T) && sq_distance <= si.switching_function.sq_cutoff)
        distance = sqrt(sq_distance)
        inv_distance = one(T) / distance

        x1 = (distance - si.R_min1) / si.σ1
        x2 = (distance - si.R_min2) / si.σ2

        # -dE/dr / r
        factor = -inv_distance^3 / T(π)^T(1.5) * (
            si.V2 * si.ΔG_free1 * exp(-x1^2) * (x1 / si.σ1 + inv_distance) / si.σ1 +
            si.V1 * si.ΔG_free2 * exp(-x2^2) * (x2 / si.σ2 + inv_distance) / si.σ2
        )

        # do we have to use the switching function (cuton <= distance <= cutoff)?
        if (sq_distance > si.switching_function.sq_cuton)
            switch_value, switch_derivative = switching_derivative(si.switching_function, sq_distance)

            # First, multiply the current force with the switching function
            factor *= switch_value

            # Second, we add the product of the energy and the derivative
            # of the switching function (the total force is the derivative of
            # a product of functions)
            factor -= switch_derivative * _eef1_energy(si, distance)
        end

        force = factor * direction

        si.a1.F += force
        si.a2.F -= force
    end
end

function compute_forces!(nbc::CharmmNonBondedComponent)
    if isempty(nbc.ff.constrained_atoms)
        map(compute_forces!, nbc.lj_interactions)
        map(compute_forces!, nbc.electrostatic_interactions)
        map(compute_forces!, nbc.solvation_interactions)
    else
        constrained_ids = getproperty.(getindex.(Ref(atoms(nbc.ff.system)), nbc.ff.constrained_atoms), :idx)
        filter_pairs = s -> (p for p in s if (p.a1.idx ∉ constrained_ids) || (p.a2.idx ∉ constrained_ids))

        map(compute_forces!, filter_pairs(nbc.lj_interactions))
        map(compute_forces!, filter_pairs(nbc.electrostatic_interactions))
        map(compute_forces!, filter_pairs(nbc.solvation_interactions))
    end
    nothing
end

function count_warnings(nbc::CharmmNonBondedComponent)
    length(nbc.unassigned_lj_interactions)
end

function print_warnings(nbc::CharmmNonBondedComponent)
    for (atom_1, atom_2) in nbc.unassigned_lj_interactions
        @warn "CharmmNonBondedComponent(): cannot find vdW parameters for "             *
                "atom types $(atom_1.atom_type)-$(atom_2.atom_type) (atoms are: "       *
                "$(get_full_name(atom_1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                "$(get_full_name(atom_2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
    end
end
