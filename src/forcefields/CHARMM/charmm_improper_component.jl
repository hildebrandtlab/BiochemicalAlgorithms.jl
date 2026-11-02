export
    QuadraticImproperTorsion,
    CharmmImproperTorsionComponent

@auto_hash_equals struct QuadraticImproperTorsion{T<:Real}
    k::T
    ϕ₀::T
    a1::Atom{T}
    a2::Atom{T}
    a3::Atom{T}
    a4::Atom{T}
end

@auto_hash_equals mutable struct CharmmImproperTorsionComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    energy::Dict{String, T}

    unassigned_impropers::Vector{Tuple{Atom{T}, Atom{T}, Atom{T}, Atom{T}}}

    improper_torsions::Vector{QuadraticImproperTorsion{T}}

    function CharmmImproperTorsionComponent{T}(ff::ForceField{T}) where {T<:Real}
        new("CharmmImproperTorsion", ff, Dict{String, T}(), [])
    end
end

# looks up the improper parameters for the given types; the first atom is the central atom
#
# this mirrors BALL: within each group, the first match is used; across groups, the last
# matching group wins. Each lookup also tries wildcards for the two inner atoms.
function _charmm_improper_parameters(improper_dict::Dict{NTuple{4, String}, V}, t1, t2, t3, t4) where V
    lookup(I, J, K, L) = @coalesce(
        get(improper_dict, (I, J,   K,   L), missing),
        get(improper_dict, (I, "*", "*", L), missing)
    )

    groups = (
        ((t1, t2, t3, t4), (t1, t2, t4, t3)),
        ((t1, t3, t2, t4), (t1, t3, t4, t2)),
        ((t1, t4, t2, t3), (t1, t4, t3, t2), (t1, "*", "*", t4))
    )

    result = missing
    for group in groups
        for key in group
            p = lookup(key...)
            if !ismissing(p)
                result = p
                break
            end
        end
    end

    result
end

function setup!(ic::CharmmImproperTorsionComponent{T}) where T
    ff = ic.ff
    empty!(ic.unassigned_impropers)

    improper_torsions = Vector{QuadraticImproperTorsion{T}}()
    ic.improper_torsions = improper_torsions

    improper_section = extract_section(ff.parameters, "ImproperTorsions")

    unit_k  = get(improper_section.properties, "unit_k",     "kcal/mol")
    unit_ϕ₀ = get(improper_section.properties, "unit_phase", "°")

    # ball used only ascii for the units; this clashes with the convention in Unitful.jl
    if unit_ϕ₀ ∈ ("deg", "degree")
        unit_ϕ₀ = "°"
    end

    k_factor  = T(ustrip((1uparse(unit_k))  |> u"kJ/mol"))
    ϕ₀_factor = T(ustrip((1uparse(unit_ϕ₀)) |> u"rad"))

    improper_dict = Dict{NTuple{4, String}, @NamedTuple{k::T, ϕ₀::T}}(
        (String(r.I), String(r.J), String(r.K), String(r.L)) => (k = k_factor * T(r.k), ϕ₀ = ϕ₀_factor * T(r.phase))
        for r in eachrow(improper_section.data)
    )

    # the atoms forming the impropers are listed per residue
    residue_impropers = haskey(ff.parameters.sections, "ResidueImproperTorsions") ?
        extract_section(ff.parameters, "ResidueImproperTorsions").data :
        DataFrame()

    if nrow(residue_impropers) == 0 || !all(c -> c ∈ names(residue_impropers), ["name", "A", "B", "C", "D"])
        push!(ff.warnings,
            "CharmmImproperTorsionComponent(): no usable section ResidueImproperTorsions in " *
            "$(ff.parameters.filename); improper torsions are disabled"
        )
        return
    end

    impropers_by_residue = Dict{String, Vector{NTuple{4, String}}}()
    for r in eachrow(residue_impropers)
        push!(
            get!(impropers_by_residue, String(r.name), NTuple{4, String}[]),
            (String(r.A), String(r.B), String(r.C), String(r.D))
        )
    end

    neighbours = _fragment_neighbours(ff.system)
    atoms_by_name(f::Fragment) = Dict{String, Atom{T}}(a.name => a for a in atoms(f))

    for frag in fragments(ff.system)
        res_name = get_full_name(frag)
        haskey(impropers_by_residue, res_name) || continue

        prev_idx, next_idx = get(neighbours, frag.idx, (nothing, nothing))

        candidates = Dict{String, Atom{T}}()
        merge!(candidates, atoms_by_name(frag))
        isnothing(prev_idx) || merge!(candidates,
            Dict("-" * n => a for (n, a) in atoms_by_name(fragment_by_idx(parent(frag), prev_idx))))
        isnothing(next_idx) || merge!(candidates,
            Dict("+" * n => a for (n, a) in atoms_by_name(fragment_by_idx(parent(frag), next_idx))))

        for names in impropers_by_residue[res_name]
            atoms_found = map(n -> get(candidates, n, nothing), names)

            if any(isnothing, atoms_found)
                push!(ff.warnings,
                    "CharmmImproperTorsionComponent(): could not find improper torsion for " *
                    "$(res_name):$(join(names, "-"))"
                )
                foreach(a -> push!(ff.unassigned_atoms, a), filter(!isnothing, atoms_found))
                continue
            end

            a1, a2, a3, a4 = atoms_found

            p = _charmm_improper_parameters(improper_dict, a1.atom_type, a2.atom_type, a3.atom_type, a4.atom_type)

            if ismissing(p)
                push!(ic.unassigned_impropers, (a1, a2, a3, a4))

                push!(ff.unassigned_atoms, a1)
                push!(ff.unassigned_atoms, a2)
                push!(ff.unassigned_atoms, a3)
                push!(ff.unassigned_atoms, a4)

                if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
                    throw(TooManyErrors())
                end
            else
                push!(improper_torsions, QuadraticImproperTorsion{T}(p.k, p.ϕ₀, a1, a2, a3, a4))
            end
        end
    end
end

function update!(::CharmmImproperTorsionComponent)
    nothing
end

# the improper angle is the angle between the planes (a1, a2, a3) and (a2, a3, a4),
# where a1 is the central atom
@inline function compute_energy(it::QuadraticImproperTorsion{T})::T where T
    ba = it.a2.r .- it.a1.r
    bc = it.a2.r .- it.a3.r
    bd = it.a2.r .- it.a4.r

    bcxba = cross(bc, ba)
    bcxbd = cross(bc, bd)

    length_bcxba = norm(bcxba)
    length_bcxbd = norm(bcxbd)

    if length_bcxba == zero(T) || length_bcxbd == zero(T)
        return zero(T)
    end

    cos_ϕ = clamp(dot(bcxba, bcxbd) / (length_bcxba * length_bcxbd), -one(T), one(T))

    it.k * (acos(cos_ϕ) - it.ϕ₀)^2
end

function compute_energy!(ic::CharmmImproperTorsionComponent{T})::T where T
    improper_torsion_energy = mapreduce(compute_energy, +, ic.improper_torsions; init=zero(T))

    ic.energy["Improper Torsion"] = improper_torsion_energy

    improper_torsion_energy
end

function compute_forces!(it::QuadraticImproperTorsion{T}) where T
    ba = it.a2.r .- it.a1.r
    bc = it.a2.r .- it.a3.r
    bd = it.a2.r .- it.a4.r

    bcxba = cross(bc, ba)
    bcxbd = cross(bc, bd)

    length_bcxba = norm(bcxba)
    length_bcxbd = norm(bcxbd)

    if length_bcxba == zero(T) || length_bcxbd == zero(T)
        return
    end

    inv_length_bcxba = one(T) / length_bcxba
    inv_length_bcxbd = one(T) / length_bcxbd

    cos_ϕ = clamp(dot(bcxba, bcxbd) * inv_length_bcxba * inv_length_bcxbd, -one(T), one(T))
    sin_ϕ = sqrt(one(T) - cos_ϕ^2)

    if sin_ϕ <= zero(T)
        return
    end

    # dE/dcos(ϕ)
    factor = -T(2) * it.k * (acos(cos_ϕ) - it.ϕ₀) / sin_ϕ

    inv_length_bcxba_2 = inv_length_bcxba^2
    inv_length_bcxbd_2 = inv_length_bcxbd^2

    bcba = dot(bc, ba)
    bcbd = dot(bc, bd)
    babd = dot(ba, bd)
    ba2  = dot(ba, ba)
    bc2  = dot(bc, bc)
    bd2  = dot(bd, bd)

    inv_lengths = inv_length_bcxba * inv_length_bcxbd

    ∂cosϕ∂bc = cos_ϕ * (inv_length_bcxbd_2 * (bcbd * bd - bd2 * bc) + inv_length_bcxba_2 * (bcba * ba - ba2 * bc)) +
        inv_lengths * (2 * babd * bc - bcba * bd - bcbd * ba)
    ∂cosϕ∂bd = cos_ϕ * inv_length_bcxbd_2 * (bcbd * bc - bc2 * bd) + inv_lengths * (bc2 * ba - bcba * bc)
    ∂cosϕ∂ba = cos_ϕ * inv_length_bcxba_2 * (bcba * bc - bc2 * ba) + inv_lengths * (bc2 * bd - bcbd * bc)

    it.a1.F += factor * ∂cosϕ∂ba
    it.a2.F -= factor * (∂cosϕ∂bc + ∂cosϕ∂bd + ∂cosϕ∂ba)
    it.a3.F += factor * ∂cosϕ∂bc
    it.a4.F += factor * ∂cosϕ∂bd
end

function compute_forces!(ic::CharmmImproperTorsionComponent)
    map(compute_forces!, ic.improper_torsions)

    nothing
end

function count_warnings(ic::CharmmImproperTorsionComponent)
    length(ic.unassigned_impropers)
end

function print_warnings(ic::CharmmImproperTorsionComponent)
    for ui in ic.unassigned_impropers
        a1, a2, a3, a4 = ui

        @warn "CharmmImproperTorsionComponent(): cannot find improper torsion parameters for "      *
            "atom types $(a1.atom_type)-$(a2.atom_type)-$(a3.atom_type)-$(a4.atom_type) (atoms are: " *
            "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a4, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
    end
end
