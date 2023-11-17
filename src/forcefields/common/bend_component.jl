export QuadraticAngleBend, QuadraticBendComponent

@auto_hash_equals struct QuadraticAngleBend{T<:Real}
    θ₀::T
    k::T
    a1::Atom{T}
    a2::Atom{T}
    a3::Atom{T}
end

@auto_hash_equals mutable struct QuadraticBendComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    unassigned_bends::Vector{Tuple{Atom{T}, Atom{T}, Atom{T}}}
    bends::Vector{QuadraticAngleBend{T}}

    function QuadraticBendComponent{T}(ff::ForceField{T}) where {T<:Real}
        new("QuadraticAngleBend", ff, Dict{Symbol, Any}(), Dict{String, T}(), [])
    end
end

function setup!(qbc::QuadraticBendComponent{T}) where {T<:Real}
    ff = qbc.ff
    
    # extract the parameter section for quadratic angle bends
    bend_section = extract_section(qbc.ff.parameters, "QuadraticAngleBend")
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

    bends_dict = Dict(
        Tuple(k) => bend_combinations[k] 
        for k in keys(bend_combinations)
    )

    bends = Vector{QuadraticAngleBend}()

    for atom in eachatom(ff.system)
        bs = bonds(atom)

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
                    get(bends_dict, (type_a1, type_a2, type_a3,), missing),
                    get(bends_dict, (type_a3, type_a2, type_a1,), missing),
                    get(bends_dict, ("*",     type_a2, "*"    ,), missing)
                )

                if ismissing(qab)
                    push!(qbc.unassigned_bends, (a1, a2, a3))
                          
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

    qbc.bends = bends
end

function update!(qbc::QuadraticBendComponent{T}) where {T<:Real}
    nothing
end

@inline function compute_energy(qab::QuadraticAngleBend{T})::T where {T<:Real}
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

function compute_energy(qbc::QuadraticBendComponent{T})::T where {T<:Real}
    # iterate over all bends in the system
    total_energy = mapreduce(compute_energy, +, qbc.bends; init=zero(T))

    qbc.energy["Angle Bends"] = total_energy

    total_energy
end

function compute_forces(qab::QuadraticAngleBend{T}) where {T<:Real}
    # calculate the vectors between the atoms and normalize if possible
    v1 = qab.a1.r - qab.a2.r
    v2 = qab.a3.r - qab.a2.r

    v1_length = norm(v1)
    v2_length = norm(v2)

    if v1_length == zero(T) || v2_length == zero(T)
        return
    end

    v1 /= v1_length
    v2 /= v2_length

    cos_θ = dot(v1, v2)
    θ = if cos_θ > T(1.0)
            zero(T)
        elseif cos_θ < T(-1.0)
            π
        else
            acos(cos_θ)
        end

    factor = T(2*force_prefactor) * qab.k * (θ - qab.θ₀)

    # calculate the cross product of v1 and v2 and normalize if possible
    crossv1v2 = normalize(cross(v1, v2))

    if isnan(crossv1v2[1])
        return
    end

    n1 = (cross(v1, crossv1v2)) * factor / v1_length
    n2 = (cross(v2, crossv1v2)) * factor / v2_length

    #@info "$(get_full_name(qab.a1))<->$(get_full_name(qab.a2))<->" *
    #      "$(get_full_name(qab.a3)) $(n1) $(n2)"

    qab.a1.F -= n1
    qab.a2.F += n1
    qab.a2.F -= n2
    qab.a3.F += n2
end

function compute_forces(qbc::QuadraticBendComponent{T}) where {T<:Real}
    # iterate over all bends in the system
    map(compute_forces, qbc.bends)

    nothing
end

function count_warnings(qbc::QuadraticBendComponent{T}) where {T<:Real}
    length(qbc.unassigned_bends)
end

function print_warnings(qbc::QuadraticBendComponent{T}) where {T<:Real}
    for ub in qbc.unassigned_bends
        a1, a2, a3 = ub

        type_a1 = a1.atom_type
        type_a2 = a2.atom_type
        type_a3 = a3.atom_type

        @warn "QuadraticBendComponent(): cannot find bend parameters for "        *
                          "atom types $(type_a1)-$(type_a2)-$(type_a3) (atoms are: "          *
                          "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                          "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                          "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
    end
end