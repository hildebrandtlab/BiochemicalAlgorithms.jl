export QuadraticBondStretch, QuadraticStretchComponent

@auto_hash_equals struct QuadraticBondStretch{T<:Real}
    r0::T
    k::T
    a1::Atom{T}
    a2::Atom{T}
end

@auto_hash_equals struct QuadraticStretchComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    stretches::AbstractVector{QuadraticBondStretch{T}}

    function QuadraticStretchComponent{T}(ff::ForceField{T}) where {T<:Real}
        # exctract the parameter section for quadratic bond stretches
        stretch_section = extract_section(ff.parameters, "QuadraticBondStretch")
        stretch_df = stretch_section.data

        unit_k  = get(stretch_section.properties, "unit_k", "kcal/mol")
        unit_r0 = get(stretch_section.properties, "unit_r0", "angstrom")

        # ball used to write Angstrom with a capitcal letter; this clashes with the convention in Unitful.jl
        if unit_r0 == "Angstrom"
            unit_r0 = "angstrom"
        end

        k_factor  = ustrip((1uparse(unit_k))  |> u"kJ/mol")
        r0_factor = ustrip((1uparse(unit_r0)) |> u"angstrom")

        # group the stretch parameters by type_i, type_j combinations
        stretch_combinations = groupby(stretch_df, ["I", "J"])

        stretches = Vector{QuadraticBondStretch}(undef, nbonds(ff.system))

        # iterate over all bonds in the system
        for (i, bond) in enumerate(bonds(ff.system))
            a1 = atom_by_idx(parent_system(ff.system), bond.a1)
            a2 = atom_by_idx(parent_system(ff.system), bond.a2)

            type_a1 = a1.atom_type
            type_a2 = a2.atom_type

            if has_flag(bond, :TYPE__HYDROGEN)
                # skip hydrogen bonds
                stretches[i] = QuadraticBondStretch(one(T), zero(T), a1, a2)
            end
       
            qbs = coalesce(
                get(stretch_combinations, (I=type_a1, J=type_a2,), missing),
                get(stretch_combinations, (I=type_a2, J=type_a1,), missing),
                get(stretch_combinations, (I=type_a1, J="*",    ), missing),
                get(stretch_combinations, (I="*",     J=type_a2,), missing),
                get(stretch_combinations, (I="*",     J="*",    ), missing)
            )

            stretches[i] = if ismissing(qbs)

                @warn "QuadraticStretchComponent(): cannot find stretch parameters for " *
                    "atom types $(type_a1)-$(type_a2) (atoms are: "                      *
                    "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"  *
                    "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"

                push!(ff.unassigned_atoms, a1)
                push!(ff.unassigned_atoms, a2)

                if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
                    throw(TooManyErrors())
                end

                # we don't want to get any force or energy component from this stretch
                QuadraticBondStretch(one(T), zero(T), a1, a2)
            else
                QuadraticBondStretch(
                    T(r0_factor*only(qbs.r0)), 
                    T(k_factor*only(qbs.k)),
                    a1,
                    a2
                )
            end
        end
        
        new("QuadraticBondStretch", ff, stretches)
    end
end

@inline function compute_energy(qbs::QuadraticBondStretch{T}) where {T<:Real}
    d = distance(qbs.a1.r, qbs.a2.r)

    qbs.k * (d - qbs.r0)^2
end

function compute_energy(qsc::QuadraticStretchComponent{T}) where {T<:Real}
    # iterate over all bonds in the system

    mapreduce(compute_energy, +, qsc.stretches; init=zero(T))
end