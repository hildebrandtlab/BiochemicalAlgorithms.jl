export
    normalize_names!

function match_name_(residue_name::AbstractString, atom_name::AbstractString, mapping::Dict{String, String})
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
    else
        # second, try wildcard match for residue names
        if haskey(mapping, "*:$(atom_name)")
            mapped_name = split(mapping["*:$(atom_name)"], ":")

            residue_name = "*"
            atom_name    = mapped_name[2]

            hit = true
        end
    end

    return hit, residue_name, atom_name
end

function do_match_(residue_name::AbstractString, residue_suffix::AbstractString, atom_name::AbstractString, mapping::Dict{String, String})
    # first, try to match exactly
    if residue_suffix != ""
        full_residue_name = residue_name * residue_suffix

        # try to match with the full residue name
        hit, res_name, atom_name = match_name_(full_residue_name, atom_name, mapping)

        if hit
            return hit, res_name, atom_name
        end
    end

    # this did not work; try to match with non-terminal residue names
    hit, res_name, atom_name = match_name_(residue_name, atom_name, mapping)

    if hit
        return hit, res_name, atom_name
    end

    # ok; let's try with a wildcard for the residue name instead
    return match_name_("*$(residue_suffix)", atom_name, mapping)
end

function get_suffix_(frag::Fragment)
    #if is_c_terminal(frag)
    if has_flag(frag, :C_TERMINAL)
        return "-C"
    #elseif is_n_terminal(frag)
    elseif has_flag(frag, :N_TERMINAL)
        return "-N"
    end

    ""
end

function count_hits_(scheme::DBNameMapping, frag::Fragment{T}) where {T<:Real}
    res_name = frag.name

    hits = 0
    res_name_suffix = get_suffix_(frag)

    for atom in atoms(frag)
        hit, res_name, _ = do_match_(res_name, res_name_suffix, atom.name, scheme.mappings)
        if hit
            hits += 1
        end
    end

    hits
end

function count_hits_(scheme::DBNameMapping, frags::FragmentTable{T}) where {T<:Real}
    if length(frags) == 0
        return
    end

    hits = 0

    for frag in frags
        hits += count_hits_(scheme, frag)
    end

    hits
end

function normalize_fragments_!(frags::FragmentTable{T}, mapping::Dict{String, String}) where {T<:Real}
    for frag in frags
        res_name = frag.name
        res_name_suffix = get_suffix_(frag)

        for atom in atoms(frag)
            hit, new_res_name, atom_name = do_match_(res_name, res_name_suffix, atom.name, mapping)

            if hit
                atom.name = atom_name

                if !startswith(new_res_name, "*")
                    frag.name = new_res_name
                end
            end
        end
    end
end

function normalize_names!(
        m::AbstractAtomContainer{T},
        fdb::FragmentDB;
        naming_standard = "") where {T<:Real}

    # start by labelling all terminal fragments to speed up terminal lookups later on
    label_terminal_fragments!(m)

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
            mapping_candidates[scheme] = count_hits_(fdb.name_mappings[scheme], frags)
        end

        # find the mapping with greatest number of hits
        hits, mapping = findmax(mapping_candidates)

        # and apply it
        if hits > 0
            normalize_fragments_!(frags, fdb.name_mappings[mapping].mappings)
        else
            @warn "normalize_names could not find a suitable mapping for $(chain)!"
        end
    end
end
