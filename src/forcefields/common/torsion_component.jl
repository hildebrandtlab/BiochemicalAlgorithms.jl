export
    CosineTorsion,
    TorsionComponent

@auto_hash_equals struct CosineTorsion{T<:Real}
    V::Vector{T}
    ϕ₀::Vector{T}
    f::Vector{Int}
    div::Vector{Int}
    a1::Atom{T}
    a1r::Vector3{T}
    a2::Atom{T}
    a2r::Vector3{T}
    a3::Atom{T}
    a3r::Vector3{T}
    a4::Atom{T}
    a4r::Vector3{T}
end

@auto_hash_equals mutable struct TorsionComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    energy::Dict{String, T}

    unassigned_torsions::Vector{Tuple{Atom{T}, Atom{T}, Atom{T}, Atom{T}, Bool}}

    proper_torsions::Vector{CosineTorsion{T}}
    improper_torsions::Vector{CosineTorsion{T}}

    function TorsionComponent{T}(ff::ForceField{T}) where {T<:Real}
        new("Torsion", ff, Dict{String, Any}(), [])
    end
end

function _get_torsion_data(ff::ForceField{T}, section::String)::Tuple{T, T, GroupedDataFrame{DataFrame}} where {T<:Real}
    # extract the parameter section for quadratic bond stretches
    torsion_section = extract_section(ff.parameters, section)
    torsion_df = torsion_section.data

    unit_V  = get(torsion_section.properties, "unit_V",    "kcal/mol")
    unit_ϕ₀ = get(torsion_section.properties, "unit_phi0", "°"       )

    # ball used only ascii for the units; this clashes with the convention in Unitful.jl
    if unit_ϕ₀ == "degree"
        unit_ϕ₀ = "°"
    end

    V_factor  = T(ustrip((1uparse(unit_V))  |> u"kJ/mol"))
    ϕ₀_factor = T(ustrip((1uparse(unit_ϕ₀)) |> u"rad"))

    # group the torsion parameters by type_i, type_j, type_k, type_l combinations
    torsion_combinations = groupby(torsion_df, ["I", "J", "K", "L"])

    V_factor, ϕ₀_factor, torsion_combinations
end

function _try_assign_torsion!(
        tc::TorsionComponent{T},
        torsions::Vector{CosineTorsion{T}},
        torsion_combinations::Dict{K, V},
        a1::Atom{T},
        a2::Atom{T},
        a3::Atom{T},
        a4::Atom{T},
        V_factor::T,
        ϕ₀_factor::T,
        is_proper::Bool) where {T<:Real, K, V}

    ff = tc.ff

    type_a1::String = a1.atom_type
    type_a2::String = a2.atom_type
    type_a3::String = a3.atom_type
    type_a4::String = a4.atom_type

    # now, search torsion parameters for (a1, a2, a3, a4)
    pt = coalesce(
        get(torsion_combinations, (type_a1, type_a2, type_a3, type_a4,), missing),
        get(torsion_combinations, (type_a4, type_a3, type_a2, type_a1,), missing),
        get(torsion_combinations, ("*",     type_a2, type_a3, type_a4,), missing),
        get(torsion_combinations, ("*",     type_a3, type_a2, type_a1,), missing),
        get(torsion_combinations, (type_a1, type_a2, type_a3, "*",    ), missing),
        get(torsion_combinations, (type_a4, type_a3, type_a2, "*",    ), missing),
        get(torsion_combinations, ("*",     type_a2, type_a3, "*"    ,), missing),
        get(torsion_combinations, ("*",     type_a3, type_a2, "*"    ,), missing),
        get(torsion_combinations, ("*",     "*",     type_a3, type_a4,), missing)
    )

    if ismissing(pt)
        push!(tc.unassigned_torsions, (a1, a2, a3, a4, is_proper))

        push!(ff.unassigned_atoms, a1)
        push!(ff.unassigned_atoms, a2)
        push!(ff.unassigned_atoms, a3)
        push!(ff.unassigned_atoms, a4)

        if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
            throw(TooManyErrors())
        end
    else
        # extract the line containing the number of frequencies in the Fourier series
        n = only(filter(d -> d.n == "N", pt)).div

        # and the lines containing the individual components
        pts = filter(d -> d.n != "N", pt)

        if n != nrow(pts)
            @warn "Invalid torsion configuration for atom types $(type_a1)-$(type_a2)-$(type_a3)-$(type_a4)!"
        else
            push!(torsions,
                CosineTorsion(
                    V_factor  .* getproperty.(pts, :V),
                    ϕ₀_factor .* getproperty.(pts, :phi0),
                    getproperty.(pts, :f),
                    getproperty.(pts, :div),
                    a1, a1.r,
                    a2, a2.r,
                    a3, a3.r,
                    a4, a4.r
                )
            )
        end
    end
end

function setup!(tc::TorsionComponent{T}) where {T<:Real}
    # first, set up the proper torsions
    _setup_proper_torsions!(tc)

    # now, set up the improper torsions
    _setup_improper_torsions!(tc)
end

function _setup_improper_torsions!(tc::TorsionComponent{T}) where {T<:Real}
    ff = tc.ff

    V_factor, ϕ₀_factor, torsion_combinations = _get_torsion_data(ff, "ImproperTorsions")

    # extract the parameter sections containing all possible improper torsion atoms
    impropers = extract_section(ff.parameters, "ResidueImproperTorsions").data

    _parse_combi_tuple(t) = (n = t.n, div = t.div, V = t.V, phi0 = t.phi0, f = t.f)
    torsion_dict = Dict(
        Tuple(k) => [_parse_combi_tuple(r) for r in eachrow(torsion_combinations[k])]
        for k in keys(torsion_combinations)
    )

    improper_torsions = Vector{CosineTorsion{T}}()

    # check for each potential improper torsion atom (every atom having three bonds)
    # whether it is contained in the list of impropers
    for atom in atoms(ff.system)
        bs = bonds(atom)
        if length(bs) == 3
            if get_full_name(atom) ∈ impropers.name
                for i_1 in eachindex(bs)
                    bond_1 = bs[i_1]
                    a3 = atom
                    a4 = get_partner(bond_1, atom)

                    for i_2 = (i_1+1):length(bs)
                        bond_2 = bs[i_2]
                        a2 = get_partner(bond_2, atom)

                        for i_3 = (i_2+1):length(bs)
                            bond_3 = bs[i_3]
                            a1 = get_partner(bond_3, atom)

                            _try_assign_torsion!(
                                tc,
                                improper_torsions,
                                torsion_dict,
                                a1, a2, a3, a4,
                                V_factor,
                                ϕ₀_factor,
                                false
                            )
                        end
                    end
                end
            end
        end
    end

    tc.improper_torsions = improper_torsions
end

function _setup_proper_torsions!(tc::TorsionComponent{T}) where {T<:Real}
    ff = tc.ff

    V_factor, ϕ₀_factor, torsion_combinations = _get_torsion_data(ff, "Torsions")

    _parse_combi_tuple(t) = (n = t.n, div = t.div, V = t.V, phi0 = t.phi0, f = t.f)
    torsion_dict = Dict{NTuple{4, String}, Vector{@NamedTuple{n::String, div::Int, V::T, phi0::T, f::Int}}}(
        Tuple(k) => [_parse_combi_tuple(r) for r in eachrow(torsion_combinations[k])]
        for k in keys(torsion_combinations)
    )

    proper_torsions = Vector{CosineTorsion{T}}()

    for atom::Atom{T} in atoms(ff.system)
        bs = bonds(atom)

        for bond_1 in bs
            if has_flag(bond_1, :TYPE__HYDROGEN)
                continue
            end

            if atom.idx == bond_1.a1
                # central atoms
                a2 = atom
                a3::Atom{T} = atom_by_idx(parent_system(atom), bond_1.a2)

                for bond_2 in bs
                    if has_flag(bond_2, :TYPE__HYDROGEN)
                        continue
                    end

                    if bond_2.a2 == bond_1.a2
                        continue
                    end

                    # determine the first atom
                    a1::Atom{T} = ((bond_2.a1 == atom.idx) ? atom_by_idx(parent_system(atom), bond_2.a2)
                                                           : atom_by_idx(parent_system(atom), bond_2.a1))

                    for bond_3 in bonds(atom_by_idx(parent_system(atom), bond_1.a2))
                        if has_flag(bond_3, :TYPE__HYDROGEN)
                            continue
                        end

                        if bond_3.a1 == a2.idx
                            continue
                        end

                        # determine the fourth atom a4
                        a4::Atom{T} = ((bond_3.a1 == a3.idx) ? atom_by_idx(parent_system(atom), bond_3.a2)
                                                             : atom_by_idx(parent_system(atom), bond_3.a1))

                        _try_assign_torsion!(
                            tc,
                            proper_torsions,
                            torsion_dict,
                            a1, a2, a3, a4,
                            V_factor,
                            ϕ₀_factor,
                            true
                        )
                    end
                end
            end
        end
    end

    tc.proper_torsions = proper_torsions
end

function update!(::TorsionComponent{T}) where {T<:Real}
    nothing
end

@inline function compute_energy(pt::CosineTorsion{T})::T where {T<:Real}
    energy = zero(T)

    a23 = pt.a3r .- pt.a2r

    cross2321 = normalize(cross(a23, pt.a1r .- pt.a2r))
    cross2334 = normalize(cross(a23, pt.a4r .- pt.a3r))

    if !isnan(cross2321[1]) && !isnan(cross2334[1])
        cos_ϕ = clamp(dot(cross2321, cross2334), T(-1.0), T(1.0))

        terms = pt.V ./ pt.div .* (Ref(1) .+ cos.(pt.f .* Ref(acos(cos_ϕ)) .- pt.ϕ₀))

        @debug "$(get_full_name(pt.a1))-$(get_full_name(pt.a2))-" *
               "$(get_full_name(pt.a3))-$(get_full_name(pt.a4)) " *
               "$(cos_ϕ) terms: $(terms)"

        energy = sum(terms)
    end

    energy
end

function compute_energy!(tc::TorsionComponent{T})::T where {T<:Real}
    # iterate over all proper torsions in the system
    proper_torsion_energy   = mapreduce(compute_energy, +, tc.proper_torsions; init=zero(T))
    improper_torsion_energy = mapreduce(compute_energy, +, tc.improper_torsions; init=zero(T))

    total_energy = proper_torsion_energy + improper_torsion_energy

    tc.energy["Proper Torsion"]   = proper_torsion_energy
    tc.energy["Improper Torsion"] = improper_torsion_energy

    total_energy
end

function compute_forces!(ct::CosineTorsion{T}) where {T<:Real}
    a21 = ct.a1r .- ct.a2r
    a23 = ct.a3r .- ct.a2r
    a34 = ct.a4r .- ct.a3r

    cross2321 = cross(a23, a21)
    cross2334 = cross(a23, a34)

    length_cross2321 = norm(cross2321)
    length_cross2334 = norm(cross2334)

    if length_cross2321 ≠ zero(T) && length_cross2334 ≠ zero(T)
        cos_ϕ = clamp(dot(cross2321, cross2334) / (length_cross2321 * length_cross2334), T(-1.0), T(1.0))

        terms = -ct.V./ct.div .* ct.f .* (sin.(ct.f .* Ref(acos(cos_ϕ)) .- ct.ϕ₀))

        # multiply with the barrier height and the factor for unit conversion
        # from kJ/(mol A) -> J/(mol m) -> N
        ∂E∂ϕ = sum(T(force_prefactor) .* terms)

        @debug "$(get_full_name(ct.a1))-$(get_full_name(ct.a2))-" *
               "$(get_full_name(ct.a3))-$(get_full_name(ct.a4)) " *
               "$(cos_ϕ) terms: $(terms) $(∂E∂ϕ)"

        direction = dot(cross(cross2321, cross2334), a23)

        if direction > 0.0
            ∂E∂ϕ *= -1
        end

        a13 = ct.a3r .- ct.a1r
        a24 = ct.a4r .- ct.a2r

        dEdt =  (∂E∂ϕ / (length_cross2321^2 * norm(a23)) * cross(cross2321, a23))
        dEdu = -(∂E∂ϕ / (length_cross2334^2 * norm(a23)) * cross(cross2334, a23))

        @debug "$(get_full_name(ct.a1))<->$(get_full_name(ct.a2))<->" *
              "$(get_full_name(ct.a3))<->$(get_full_name(ct.a4)) "   *
              "$(cross(dEdt, a23)); $(cross(a13, dEdt) + cross(dEdu, a34));" *
              "$(cross(a21, dEdt) + cross(a24, dEdu)); $(cross(dEdu, a23))"

        ct.a1.F += cross(dEdt, a23)
        ct.a2.F += cross(a13, dEdt) + cross(dEdu, a34)
        ct.a3.F += cross(a21, dEdt) + cross(a24, dEdu)
        ct.a4.F += cross(dEdu, a23)
    end
end

function compute_forces!(tc::TorsionComponent{T}) where {T<:Real}
    map(compute_forces!, tc.proper_torsions)
    map(compute_forces!, tc.improper_torsions)

    nothing
end

function count_warnings(tc::TorsionComponent{T}) where {T<:Real}
    length(tc.unassigned_torsions)
end

function print_warnings(tc::TorsionComponent{T}) where {T<:Real}
    for ut in tc.unassigned_torsions
        a1, a2, a3, a4, is_proper = ut

        type_a1::String = a1.atom_type
        type_a2::String = a2.atom_type
        type_a3::String = a3.atom_type
        type_a4::String = a4.atom_type

        torsion_type = (is_proper) ? "proper" : "improper"

        @warn "TorsionComponent(): cannot find $(torsion_type) torsion parameters for "      *
            "atom types $(type_a1)-$(type_a2)-$(type_a3)-$(type_a4) (atoms are: " *
            "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a4, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
    end
end
