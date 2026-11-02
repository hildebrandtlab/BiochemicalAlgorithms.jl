export
    CharmmTorsionComponent

@auto_hash_equals mutable struct CharmmTorsionComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    energy::Dict{String, T}

    unassigned_torsions::Vector{Tuple{Atom{T}, Atom{T}, Atom{T}, Atom{T}, Bool}}

    proper_torsions::Vector{CosineTorsion{T}}

    function CharmmTorsionComponent{T}(ff::ForceField{T}) where {T<:Real}
        new("CharmmTorsion", ff, Dict{String, T}(), [])
    end
end

# CHARMM only uses exact matches and wildcards for both outer atoms
_charmm_torsion_candidate_keys(t1::String, t2::String, t3::String, t4::String) = (
    (t1,  t2, t3, t4 ),
    (t4,  t3, t2, t1 ),
    ("*", t2, t3, "*"),
    ("*", t3, t2, "*")
)

# maps each fragment idx to the idx of the previous and next fragment in its chain
# (or `nothing` if there is no such fragment)
function _fragment_neighbours(ac::AbstractAtomContainer)
    neighbours = Dict{Int, Tuple{MaybeInt, MaybeInt}}()

    for c in chains(ac)
        fs = fragments(c)
        for i in eachindex(fs)
            neighbours[fs[i].idx] = (
                i > firstindex(fs) ? fs[i-1].idx : nothing,
                i < lastindex(fs)  ? fs[i+1].idx : nothing
            )
        end
    end

    neighbours
end

# returns the name of `atom` as seen from fragment `ref`, i.e., prefixed with "-" ("+") if
# the atom belongs to the previous (next) fragment in the chain, or with "=" if it belongs
# to a different fragment connected to `ref` via a disulphide bond
function _relative_atom_name(atom::Atom, ref::Fragment, neighbours::Dict{Int, Tuple{MaybeInt, MaybeInt}})
    frag = parent_fragment(atom)

    if isnothing(frag) || frag.idx == ref.idx
        return atom.name
    end

    prev_idx, next_idx = get(neighbours, ref.idx, (nothing, nothing))

    if frag.idx == prev_idx
        return "-" * atom.name
    end

    if frag.idx == next_idx
        return "+" * atom.name
    end

    if _has_ssbond(ref) && _has_ssbond(frag)
        return "=" * atom.name
    end

    atom.name
end

function setup!(tc::CharmmTorsionComponent{T}) where T
    ff = tc.ff
    empty!(tc.unassigned_torsions)

    V_factor, ϕ₀_factor, torsion_combinations = _get_torsion_data(ff, "Torsions")

    _parse_combi_tuple(t) = (n = t.n, div = t.div, V = t.V, phi0 = t.phi0, f = t.f)
    torsion_dict = Dict{NTuple{4, String}, Vector{@NamedTuple{n::String, div::Int, V::T, phi0::T, f::Int}}}(
        Tuple(k) => [_parse_combi_tuple(r) for r in eachrow(torsion_combinations[k])]
        for k in keys(torsion_combinations)
    )

    # if the parameter file lists the torsions for each residue, we only consider those;
    # otherwise, all torsions are computed from the connectivity
    residue_torsions = Set{NTuple{5, String}}()
    if haskey(ff.parameters.sections, "ResidueTorsions")
        for r in eachrow(extract_section(ff.parameters, "ResidueTorsions").data)
            push!(residue_torsions, (String(r.name), String(r.A), String(r.B), String(r.C), String(r.D)))
        end
    end

    neighbours = _fragment_neighbours(ff.system)

    function _is_listed(a1, a2, a3, a4)
        frag = parent_fragment(a2)

        if isnothing(frag)
            push!(ff.warnings,
                "CharmmTorsionComponent(): cannot identify torsion for atoms without fragment " *
                "$(a1.name)/$(a2.name)/$(a3.name)/$(a4.name)"
            )
            return false
        end

        res_name = get_full_name(frag)
        n1 = _relative_atom_name(a1, frag, neighbours)
        n3 = _relative_atom_name(a3, frag, neighbours)
        n4 = _relative_atom_name(a4, frag, neighbours)

        (res_name, n1, a2.name, n3, n4) ∈ residue_torsions ||
            (res_name, n4, n3, a2.name, n1) ∈ residue_torsions
    end

    proper_torsions = Vector{CosineTorsion{T}}()

    _foreach_proper_torsion(ff.system) do a1, a2, a3, a4
        if isempty(residue_torsions) || _is_listed(a1, a2, a3, a4)
            _try_assign_torsion!(
                tc,
                proper_torsions,
                torsion_dict,
                a1, a2, a3, a4,
                V_factor,
                ϕ₀_factor,
                true;
                candidate_keys = _charmm_torsion_candidate_keys
            )
        end
    end

    tc.proper_torsions = proper_torsions
end

function update!(::CharmmTorsionComponent)
    nothing
end

function compute_energy!(tc::CharmmTorsionComponent{T})::T where T
    proper_torsion_energy = mapreduce(compute_energy, +, tc.proper_torsions; init=zero(T))

    tc.energy["Proper Torsion"] = proper_torsion_energy

    proper_torsion_energy
end

function compute_forces!(tc::CharmmTorsionComponent)
    map(compute_forces!, tc.proper_torsions)

    nothing
end

function count_warnings(tc::CharmmTorsionComponent)
    length(tc.unassigned_torsions)
end

function print_warnings(tc::CharmmTorsionComponent)
    for ut in tc.unassigned_torsions
        a1, a2, a3, a4, _ = ut

        @warn "CharmmTorsionComponent(): cannot find torsion parameters for "                   *
            "atom types $(a1.atom_type)-$(a2.atom_type)-$(a3.atom_type)-$(a4.atom_type) (atoms are: " *
            "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a3, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"   *
            "$(get_full_name(a4, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"
    end
end
