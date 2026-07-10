export
    vdw_radius,
    assign_radii!,
    assign_ball_radii!,
    load_ball_radii_table

const _DEFAULT_VDW_RADIUS = 1.5

"""
    $(TYPEDSIGNATURES)

van-der-Waals radius (in Å) for the given element, as type `T`.

Values are taken from Mendeleev's element data (Bondi-style). For elements
without a published value (e.g., very heavy synthetic atoms) or for
`Elements.Unknown`, a fallback of $(_DEFAULT_VDW_RADIUS) Å is returned.
"""
@inline function vdw_radius(::Type{T}, e::ElementType) where {T<:Real}
    sym = Symbol(e)
    haskey(chem_elements, sym) || return T(_DEFAULT_VDW_RADIUS)
    val = chem_elements[sym].vdw_radius
    ismissing(val) ? T(_DEFAULT_VDW_RADIUS) : T(ustrip(u"Å", val))
end

"""
    $(TYPEDSIGNATURES)

Assign van-der-Waals radii to the atoms of `ac`. Atoms whose `radius` field is
already non-zero are left untouched unless `overwrite=true`. Returns `ac`.
"""
function assign_radii!(ac::AbstractAtomContainer{T}; overwrite::Bool=false) where T
    at = atoms(ac)
    @inbounds for i in eachindex(at.radius)
        if overwrite || iszero(at.radius[i])
            at.radius[i] = vdw_radius(T, at.element[i])
        end
    end
    ac
end

"""
    load_ball_radii_table(path::AbstractString) -> Dict{Tuple{String,String}, Float64}

Load BALL's `.siz` radii table format (e.g., `PARSE.siz`, `amber94.siz`).
Each entry is `RESIDUE:ATOM_NAME RADIUS`. Wildcard `*` matches any residue.
Returns a dictionary keyed by `(residue, atom_name)`.
"""
function load_ball_radii_table(path::AbstractString)
    rules = Dict{Tuple{String,String}, Float64}()
    for line in eachline(path)
        ln = strip(line)
        (isempty(ln) || startswith(ln, "#")) && continue
        parts = split(ln)
        length(parts) >= 2 || continue
        ra = split(parts[1], ":")
        length(ra) == 2 || continue
        rules[(String(ra[1]), String(ra[2]))] = parse(Float64, parts[2])
    end
    return rules
end

"""
    assign_ball_radii!(ac, rules::Dict{Tuple{String,String}, Float64};
                       fallback_element_radii::Bool = false)

Assign atom radii from a BALL `.siz`-style table loaded by
[`load_ball_radii_table`](@ref). Each atom is looked up by
`(residue.name, atom.name)`; if no match, falls back to wildcard
`("*", atom_name)`.

When `fallback_element_radii=false` (the default), atoms that match
neither the explicit residue rule nor the wildcard are LEFT WITH
THEIR EXISTING RADIUS (which is `0` for atoms that have never been
assigned). BALL's `radiusRuleProcessor` follows this convention —
atoms with `radius == 0` are excluded from molecular-surface
computations (see [`compute_reduced_surface`](@ref)).

Set `fallback_element_radii=true` to fall back to element-based vdW
defaults (Bondi via Mendeleev) for unmatched atoms. This INCLUDES
heteroatoms (heme Fe, exotic ions, etc.) in the surface — useful
when you want a full molecular surface, but **breaks bit-for-bit
parity with BALL** on PDBs containing HETATM records.
"""
function assign_ball_radii!(ac::AbstractAtomContainer{T},
                              rules::Dict{Tuple{String,String}, Float64};
                              fallback_element_radii::Bool = false) where T
    sentinel = -1.0  # missing-from-table marker (lets 0.0 stay distinguishable)
    for a in atoms(ac)
        pf = parent_fragment(a)
        residue = pf === nothing ? "" : string(pf.name)
        atom_name = string(a.name)
        r = get(rules, (residue, atom_name),
                get(rules, ("*", atom_name), sentinel))
        if r == sentinel
            if fallback_element_radii
                a.radius = vdw_radius(T, a.element)
            end
        else
            a.radius = T(r)
        end
    end
    ac
end
