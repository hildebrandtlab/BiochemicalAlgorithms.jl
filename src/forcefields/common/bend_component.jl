export QuadraticAngleBend, QuadraticBendComponent

@auto_hash_equals struct QuadraticAngleBend{T<:Real}
    θ₀::T
    k::T
    a1::Atom{T}
    a2::Atom{T}
    a3::Atom{T}
end

@auto_hash_equals struct QuadraticBendComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    bends::AbstractVector{QuadraticAngleBend{T}}

    function QuadraticBendComponent{T}(ff::ForceField{T}) where {T<:Real}
        # extract the parameter section for quadratic angle bends
        bend_section = extract_section(ff.parameters, "QuadraticAngleBend")
        bend_df = bend_section.data

        unit_k  = get(bend_section.properties, "unit_k", "kcal/mol")
        unit_θ₀ = get(bend_section.properties, "unit_theta0", "°")


        # ball used only ascii for the units; this clashes with the convention in Unitful.jl
        if unit_θ₀ == "degree"
            unit_θ₀ = "°"
        end

        k_factor  = ustrip((1uparse(unit_k))  |> u"kJ/mol")
        θ₀_factor = ustrip((1uparse(unit_θ₀)) |> u"rad")

        # group the bend parameters by type_i, type_j, type_k combinations
        bend_combinations = groupby(bend_df, ["I", "J", "K"])

        bends = Vector{QuadraticAngleBend}()

        for atom in atoms(ff.system)
            bs = collect(bonds(atom))

            for i in eachindex(bs)
                bond_1 = bs[i]

                if has_flag(bond_1, :TYPE__HYDROGEN)
                    continue
                end

                for j in i+1:length(bs)
                    bond_2 = bs[j]
                
                    if has_flag(bond_2, :TYPE__HYDROGEN)
                        continue
                    end

                    a1 = get_partner(bond_1, atom)
                    a2 = atom
                    a3 = get_partner(bond_2, atom)

                    type_a1 = a1.atom_type
                    type_a2 = a2.atom_type
                    type_a3 = a3.atom_type

                    qab = coalesce(
                        get(bend_combinations, (I=type_a1, J=type_a2, K=type_a3,), missing),
                        get(bend_combinations, (I=type_a3, J=type_a2, K=type_a1,), missing),
                        get(bend_combinations, (I="*",     J=type_a2, K="*"    ,), missing)
                    )

                    if ismissing(qab)
                        @warn "QuadraticBendComponent(): cannot find bend parameters for "        *
                              "atom types $(type_a1)-$(type_a2)-$(type_a3) (atoms are: "          *
                              "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                              "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                              "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
                              
                        push!(ff.unassigned_atoms, a1)
                        push!(ff.unassigned_atoms, a2)
                        push!(ff.unassigned_atoms, a3)

                        if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
                            throw(TooManyErrors())
                        end
                    else
                        push!(bends, 
                            QuadraticAngleBend(
                                T(θ₀_factor*only(qab.theta0)),
                                T(k_factor*only(qab.k)),
                                a1,
                                a2,
                                a3
                            ))
                    end
                end
            end
        end

        new("QuadraticAngleBend", ff, bends)
    end
end


@inline function compute_energy(qab::QuadraticAngleBend{T}) where {T<:Real}
    v1 = qab.a1.r - qab.a2.r
    v2 = qab.a3.r - qab.a2.r

    sq_length = squared_norm(v1) * squared_norm(v2)

    if iszero(sq_length)
        return zero(T)
    end

    cos_θ = dot(v1, v2) / sqrt(sq_length)
    θ = if cos_θ > T(1.0)
            zero(T)
        elseif cos_θ < T(-1.0)
            π
        else
            acos(cos_θ)
        end

    qab.k * (θ - qab.θ₀)^2
end

function compute_energy(qbc::QuadraticBendComponent{T}) where {T<:Real}
    # iterate over all bends in the system
    mapreduce(compute_energy, +, qbc.bends; init=zero(T))
end