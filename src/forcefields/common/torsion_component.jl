export CosineTorsion, TorsionComponent

@auto_hash_equals struct CosineTorsion{T<:Real}
    V::AbstractVector{T}
    ϕ₀::AbstractVector{T}
    f::AbstractVector{Int}
    div::AbstractVector{Int}
    a1::Atom{T}
    a2::Atom{T}
    a3::Atom{T}
    a4::Atom{T}
end

function _get_torsion_data(ff::ForceField{T}, section::String) where {T<:Real}
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
        ff::ForceField{T}, 
        torsions::Vector{CosineTorsion},
        torsion_combinations,
        a1::Atom{T}, 
        a2::Atom{T}, 
        a3::Atom{T}, 
        a4::Atom{T},
        V_factor::T,
        ϕ₀_factor::T) where {T<:Real}
    type_a1 = a1.atom_type
    type_a2 = a2.atom_type
    type_a3 = a3.atom_type
    type_a4 = a4.atom_type

    # now, search torsion parameters for (a1, a2, a3, a4)
    pt = coalesce(
        get(torsion_combinations, (I=type_a1, J=type_a2, K=type_a3, L=type_a4,), missing),
        get(torsion_combinations, (I=type_a4, J=type_a3, K=type_a2, L=type_a1,), missing),
        get(torsion_combinations, (I="*",     J=type_a2, K=type_a3, L=type_a4,), missing),
        get(torsion_combinations, (I="*",     J=type_a3, K=type_a2, L=type_a1,), missing),
        get(torsion_combinations, (I=type_a1, J=type_a2, K=type_a3, L="*",    ), missing),
        get(torsion_combinations, (I=type_a4, J=type_a3, K=type_a2, L="*",    ), missing),
        get(torsion_combinations, (I="*",     J=type_a2, K=type_a3, L="*"    ,), missing),
        get(torsion_combinations, (I="*",     J=type_a3, K=type_a2, L="*"    ,), missing),
        get(torsion_combinations, (I="*",     J="*",     K=type_a3, L=type_a4,), missing)
    )

    if ismissing(pt)
        @warn "TorsionComponent(): cannot find torsion parameters for "      *
            "atom types $(type_a1)-$(type_a2)-$(type_a3)-$(type_a4) (atoms are: " *
            "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a4, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
            
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
                    T.(V_factor  * pts.V),
                    T.(ϕ₀_factor * pts.phi0),
                    Int.(pts.f),
                    Int.(pts.div),
                    a1,
                    a2,
                    a3,
                    a4
                )
            )
        end
    end
end

@auto_hash_equals mutable struct TorsionComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    proper_torsions::AbstractVector{CosineTorsion{T}}
    improper_torsions::AbstractVector{CosineTorsion{T}}

    function TorsionComponent{T}(ff::ForceField{T}) where {T<:Real}
        this = new("Torsion", ff, Dict{Symbol, Any}(), Dict{String, Any}())

        setup!(this)
        update!(this)

        this
    end
end

function setup!(tc::TorsionComponent{T}) where {T<:Real}
    # first, set up the proper torsions
    V_factor, ϕ₀_factor, torsion_combinations = _get_torsion_data(tc.ff, "Torsions")

    # remember those parts that stay constant when only the system is updated
    tc.cache[:proper_V_factor]  = T(V_factor)
    tc.cache[:proper_ϕ₀_factor] = T(ϕ₀_factor)

    tc.cache[:proper_torsion_combinations] = torsion_combinations

    # now, set up the improper torsions
    V_factor, ϕ₀_factor, torsion_combinations = _get_torsion_data(tc.ff, "ImproperTorsions")

    # extract the parameter sections containing all possible improper torsion atoms
    impropers = extract_section(tc.ff.parameters, "ResidueImproperTorsions").data

    # and again, remember those parts that stay constant when only the system is updated
    tc.cache[:improper_V_factor]  = T(V_factor)
    tc.cache[:improper_ϕ₀_factor] = T(ϕ₀_factor)

    tc.cache[:improper_torsion_combinations] = torsion_combinations

    tc.cache[:impropers] = impropers
end

function update!(tc::TorsionComponent{T}) where {T<:Real}
    _update_proper_torsions!(tc)
    _update_improper_torsions!(tc)
end

function _update_proper_torsions!(tc::TorsionComponent{T}) where {T<:Real}
    ff = tc.ff

    torsion_combinations = tc.cache[:proper_torsion_combinations]

    V_factor  = tc.cache[:proper_V_factor]
    ϕ₀_factor = tc.cache[:proper_ϕ₀_factor]

    proper_torsions = Vector{CosineTorsion}()

    for atom in atoms(ff.system)
        bs = bonds(atom)

        for bond_1 in bs
            if has_flag(bond_1, :TYPE__HYDROGEN)
                continue
            end

            if atom.idx == bond_1.a1
                # central atoms
                a2 = atom
                a3 = atom_by_idx(parent_system(atom), bond_1.a2)

                for bond_2 in bs
                    if has_flag(bond_2, :TYPE__HYDROGEN)
                        continue
                    end

                    if bond_2.a2 == bond_1.a2
                        continue
                    end

                    # determine the first atom
                    a1 = ((bond_2.a1 == atom.idx) ? atom_by_idx(parent_system(atom), bond_2.a2)
                                                : atom_by_idx(parent_system(atom), bond_2.a1))

                    for bond_3 in bonds(atom_by_idx(parent_system(atom), bond_1.a2))
                        if has_flag(bond_3, :TYPE__HYDROGEN)
                            continue
                        end

                        if bond_3.a1 == a2.idx
                            continue
                        end
                    
                        # determine the fourth atom a4
                        a4 = ((bond_3.a1 == a3.idx) ? atom_by_idx(parent_system(atom), bond_3.a2)
                                                    : atom_by_idx(parent_system(atom), bond_3.a1))

                        _try_assign_torsion!(ff, proper_torsions, torsion_combinations, a1, a2, a3, a4, V_factor, ϕ₀_factor)
                    end
                end
            end
        end
    end

    tc.proper_torsions = proper_torsions
end

function _update_improper_torsions!(tc::TorsionComponent{T}) where {T<:Real}
    ff = tc.ff

    torsion_combinations = tc.cache[:improper_torsion_combinations]

    V_factor  = tc.cache[:improper_V_factor]
    ϕ₀_factor = tc.cache[:improper_ϕ₀_factor]

    impropers = tc.cache[:impropers]

    improper_torsions = Vector{CosineTorsion}()

    # check for each potential improper torsion atom (every atom having three bonds)
    # whether it is contained in the list of impropers
    for atom in atoms(ff.system)
        if nbonds(atom) == 3
            if get_full_name(atom) ∈ impropers.name
                bs = collect(bonds(atom))
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

                            _try_assign_torsion!(ff, improper_torsions, torsion_combinations, a1, a2, a3, a4, V_factor, ϕ₀_factor)
                        end
                    end
                end
            end
        end
    end

    tc.improper_torsions = improper_torsions
end

@inline function compute_energy(pt::CosineTorsion{T}) where {T<:Real}
    energy = zero(T)

    a21 = pt.a1.r - pt.a2.r
    a23 = pt.a3.r - pt.a2.r
    a34 = pt.a4.r - pt.a3.r

    cross2321 = normalize(cross(a23, a21))
    cross2334 = normalize(cross(a23, a34))

    if !isnan(cross2321[1]) && !isnan(cross2334[1])
        cos_ϕ = clamp(dot(cross2321, cross2334), T(-1.0), T(1.0))

        terms = pt.V./pt.div .* (1 .+ cos.(pt.f .* Ref(acos(cos_ϕ)) .- pt.ϕ₀))

        @debug "$(get_full_name(pt.a1))-$(get_full_name(pt.a2))-" *
               "$(get_full_name(pt.a3))-$(get_full_name(pt.a4)) " *
               "$(cos_ϕ) terms: $(terms)"

        energy = sum(terms)
    end

    energy
end

function compute_energy(tc::TorsionComponent{T}) where {T<:Real}
    # iterate over all proper torsions in the system
    proper_torsion_energy   = mapreduce(compute_energy, +, tc.proper_torsions; init=zero(T))
    improper_torsion_energy = mapreduce(compute_energy, +, tc.improper_torsions; init=zero(T))

    total_energy = proper_torsion_energy + improper_torsion_energy

    tc.energy["Proper Torsion"]   = proper_torsion_energy
    tc.energy["Improper Torsion"] = improper_torsion_energy

    total_energy
end