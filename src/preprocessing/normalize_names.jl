export
    normalize_names!

function _match_name(residue_name::AbstractString, atom_name::AbstractString, mapping::Dict{String, String})
    # start with the residue name
    residue_name = strip(residue_name)

    # does this residue have a different name in the mapped naming scheme?
    if haskey(mapping, "$(residue_name):*")
        residue_name = split(mapping["$(residue_name):*"], ":")[1]
    end

    # and now the atom name
    atom_name = strip(atom_name)

    hit = false

    # first, try to match exactly
    if haskey(mapping, "$(residue_name):$(atom_name)")
        mapped_name = split(mapping["$(residue_name):$(atom_name)"], ":")

        residue_name = mapped_name[1]
        atom_name    = mapped_name[2]

        hit = true

    # second, try wildcard match for residue names
    elseif haskey(mapping, "*:$(atom_name)")
        mapped_name = split(mapping["*:$(atom_name)"], ":")

        residue_name = "*"
        atom_name    = mapped_name[2]

        hit = true
    end

    return hit, residue_name, atom_name
end

function _count_hits(scheme::DBNameMapping, frag::Fragment)
    res_name = frag.name

    hits = 0
    for atom in atoms(frag)
        hit, res_name, _ = _match_name(res_name, atom.name, scheme.mappings)
        if hit
            hits += 1
        end
    end

    hits
end

@inline function _count_hits(scheme::DBNameMapping, frags::FragmentTable)
    sum(_count_hits(scheme, frag) for frag in frags; init = 0)
end

function _normalize_fragments!(frags::FragmentTable{T}, mapping::Dict{String, String}) where {T<:Real}
    for frag in frags
        res_name = frag.name

        for atom in atoms(frag)
            hit, new_res_name, atom_name = _match_name(res_name, atom.name, mapping)

            if hit
                atom.name = atom_name

                if !startswith(new_res_name, "*")
                    frag.name = new_res_name
                end
            end
        end
    end
end

"""
    normalize_names!(::AbstractAtomContainer{Float32})
    normalize_names!(::AbstractAtomContainer{T}, ::FragmentDB{T})

Attempts to normalize fragment and atom names according to the default/given
fragment database.

# Supported keyword arguments
 - `naming_standard::AbstractString = ""`
   The naming standard to be used for normalization. If empty, the default
   naming standard of the fragment database is used.
"""
function normalize_names!(
    m::AbstractAtomContainer{T},
    fdb::FragmentDB{T};
    naming_standard::AbstractString = ""
) where T

    if naming_standard == ""
        naming_standard = fdb.defaults["Naming"]
    end

    # Determine the most likely naming standard used in the given molecule
    mapping_candidates = Dict{String, Int}(
        n => 0 for (n,v) in filter(
            ((mk, mv),) -> mv.maps_to == naming_standard,
            fdb.name_mappings)
    )

    if length(mapping_candidates) == 0
        throw(ArgumentError("normalize_names: no mapping candidates for naming standard $(naming_standard)"))
    end

    # Normalization happens on a per-fragment - basis; if this Molecule does not
    # have any fragments, we do nothing
    if nfragments(m) == 0
        return
    end

    # Now, sort the fragments into parent containers, i.e., protein chains or nucleid acids,
    # if available. Normalizing names over a whole collection of fragments should give a more
    # stable estimate of the source naming scheme.
    parent_chains = chains(m)

    for chain in parent_chains
        frags = fragments(chain)

        for scheme in keys(mapping_candidates)
            mapping_candidates[scheme] = _count_hits(fdb.name_mappings[scheme], frags)
        end

        # find the mapping with greatest number of hits
        hits, mapping = findmax(mapping_candidates)

        # and apply it
        if hits > 0
            _normalize_fragments!(frags, fdb.name_mappings[mapping].mappings)
        else
            @warn "normalize_names could not find a suitable mapping for $(chain)!"
        end
    end

    nothing
end

@inline function normalize_names!(m::AbstractAtomContainer{Float32})
    normalize_names!(m, default_fragmentdb())
end
