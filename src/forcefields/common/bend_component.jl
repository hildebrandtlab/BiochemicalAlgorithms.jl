export
    compute_bend_energy,
    compute_bend_force,
    QuadraticAngleBend,
    QuadraticBendComponent

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

    function QuadraticBendComponent{T}(ff::ForceField{T}) where T
        new("QuadraticAngleBend", ff, Dict{Symbol, Any}(), Dict{String, T}(), [])
    end
end

function setup!(qbc::QuadraticBendComponent{T}) where T
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
        NTuple{3, String}(k) => bend_combinations[k]
        for k in keys(bend_combinations)
    )

    bends = Vector{QuadraticAngleBend{T}}()

    for atom in atoms(ff.system)
        bs = non_hydrogen_bonds(atom)

        for i in eachindex(bs)
            bond_1 = bs[i]

            for j in i+1:length(bs)
                bond_2 = bs[j]

                a1 = get_partner(bond_1, atom)
                a2 = atom
                a3 = get_partner(bond_2, atom)

                type_a1 = a1.atom_type
                type_a2 = a2.atom_type
                type_a3 = a3.atom_type

                qab = @coalesce(
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
                        QuadraticAngleBend{T}(
                            T(θ₀_factor*only(qab.theta0)),
                            T(k_factor*only(qab.k)),
                            a1,
                            a2,
                            a3,
                        ))
                end
            end
        end
    end

    qbc.bends = bends
end

function update!(::QuadraticBendComponent)
    nothing
end

@inline function compute_bend_energy(k::T, θ0::T, v1::Vector3{T}, v2::Vector3{T})::T where T
    sq_length = squared_norm(v1) * squared_norm(v2)

    if iszero(sq_length)
        return zero(T)
    end

    cos_θ = dot(v1, v2) / sqrt(sq_length)
    θ = if cos_θ > one(T)
            zero(T)
        elseif cos_θ < -one(T)
            π
        else
            acos(cos_θ)
        end

    k * (θ - θ0)^2
end

@inline function compute_bend_force(k::T, θ0::T, v1::Vector3{T}, v2::Vector3{T})::Tuple{Vector3{T}, Vector3{T}, Vector3{T}} where T
    v1_length = norm(v1)
    v2_length = norm(v2)

    if v1_length == zero(T) || v2_length == zero(T)
        return zero(v1), zero(v1), zero(v1)
    end

    v1_norm = v1 / v1_length
    v2_norm = v2 / v2_length

    cos_θ = dot(v1_norm, v2_norm)
    θ = if cos_θ > one(T)
            zero(T)
        elseif cos_θ < -one(T)
            π
        else
            acos(cos_θ)
        end

    factor = T(2) * k * (θ - θ0)

    crossv1v2 = normalize(cross(v1_norm, v2_norm))

    if isnan(crossv1v2[1])
        return zero(v1), zero(v1), zero(v1)
    end

    n1 = cross(v1_norm, crossv1v2) .* factor ./ v1_length
    n2 = cross(v2_norm, crossv1v2) .* factor ./ v2_length

    f1 = -n1
    f2 = n1 - n2
    f3 = n2

    f1, f2, f3
end

@inline function compute_energy(qab::QuadraticAngleBend{T})::T where T
    v1 = qab.a1.r .- qab.a2.r
    v2 = qab.a3.r .- qab.a2.r
    _compute_bend_energy(qab.k, qab.θ₀, v1, v2)
end

function compute_energy!(qbc::QuadraticBendComponent{T})::T where T
    # iterate over all bends in the system
    total_energy = mapreduce(compute_energy, +, qbc.bends; init=zero(T))

    qbc.energy["Angle Bends"] = total_energy

    total_energy
end

function compute_forces!(qab::QuadraticAngleBend{T}) where {T<:Real}
    v1 = qab.a1.r .- qab.a2.r
    v2 = qab.a3.r .- qab.a2.r
    
    f1, f2, f3 = _compute_bend_force(qab.k, qab.θ₀, v1, v2)
    
    qab.a1.F += f1
    qab.a2.F += f2
    qab.a3.F += f3
end

function compute_forces!(qbc::QuadraticBendComponent)
    # iterate over all bends in the system
    map(compute_forces!, qbc.bends)

    nothing
end

function count_warnings(qbc::QuadraticBendComponent)
    length(qbc.unassigned_bends)
end

function print_warnings(qbc::QuadraticBendComponent)
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
