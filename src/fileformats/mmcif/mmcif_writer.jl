using Printf

using BiochemicalAlgorithms

"""
    write_mmcif_impl(io::IO, ac::AbstractAtomContainer{T})

Write an atom container as PDBx/mmCIF format to the given IO stream.
"""
function write_mmcif_impl(io::IO, ac::AbstractAtomContainer{T}) where T
    name = replace(ac.name, r"\s+" => "_")
    isempty(name) && (name = "unnamed")

    println(io, "data_", name)
    println(io, "#")

    _write_atom_site(io, ac)
    _write_struct_conn(io, ac)
    _write_struct_conf(io, ac)
    _write_struct_sheet_range(io, ac)
end

# ─── CIF value quoting ───────────────────────────────────────────────

"""Quote a string value for CIF output."""
function _cif_quote(s::AbstractString)
    s = string(s)
    isempty(s) && return "."
    # Multiline values must use semicolon text blocks
    if occursin('\n', s) || occursin('\r', s)
        return ";\n$s\n;"
    end
    # No quoting needed for simple values
    if !any(c -> isspace(c), s) && !startswith(s, '_') && !startswith(s, '#') &&
       !startswith(s, '\'') && !startswith(s, '"') && s != "." && s != "?"
        return s
    end
    # Use single quotes if possible
    if !occursin('\'', s)
        return "'$s'"
    end
    # Use double quotes
    if !occursin('"', s)
        return "\"$s\""
    end
    # Fall back to semicolon text block
    return ";\n$s\n;"
end

# ─── _atom_site ───────────────────────────────────────────────────────

function _write_atom_site(io::IO, ac::AbstractAtomContainer{T}) where T
    all_atoms = atoms(ac)
    isempty(all_atoms) && return

    println(io, "loop_")
    println(io, "_atom_site.group_PDB")
    println(io, "_atom_site.id")
    println(io, "_atom_site.type_symbol")
    println(io, "_atom_site.label_atom_id")
    println(io, "_atom_site.label_alt_id")
    println(io, "_atom_site.label_comp_id")
    println(io, "_atom_site.label_asym_id")
    println(io, "_atom_site.label_entity_id")
    println(io, "_atom_site.label_seq_id")
    println(io, "_atom_site.pdbx_PDB_ins_code")
    println(io, "_atom_site.Cartn_x")
    println(io, "_atom_site.Cartn_y")
    println(io, "_atom_site.Cartn_z")
    println(io, "_atom_site.occupancy")
    println(io, "_atom_site.B_iso_or_equiv")
    println(io, "_atom_site.pdbx_formal_charge")
    println(io, "_atom_site.auth_seq_id")
    println(io, "_atom_site.auth_comp_id")
    println(io, "_atom_site.auth_asym_id")
    println(io, "_atom_site.auth_atom_id")
    println(io, "_atom_site.pdbx_PDB_model_num")

    # Build entity ID map (one entity per unique chain name)
    entity_map = Dict{String, Int}()
    entity_count = 0
    for c in chains(ac)
        if !haskey(entity_map, c.name)
            entity_count += 1
            entity_map[c.name] = entity_count
        end
    end

    for a in all_atoms
        frag = parent_fragment(a)
        chain = isnothing(frag) ? nothing : parent_chain(frag)

        group = is_hetero_atom(a) ? "HETATM" : "ATOM"
        type_sym = string(a.element)
        atom_name = _cif_quote(a.name)
        alt_id = "."

        comp_id = isnothing(frag) ? "UNK" : frag.name
        chain_id = isnothing(chain) ? "." : chain.name
        entity_id = isnothing(chain) ? 1 : get(entity_map, chain.name, 1)
        seq_id = isnothing(frag) ? 0 : frag.number
        ins_code = isnothing(frag) ? "?" : _cif_quote(get_property(frag, :insertion_code, "?"))

        x = @sprintf("%.3f", a.r[1])
        y = @sprintf("%.3f", a.r[2])
        z = @sprintf("%.3f", a.r[3])

        occ = @sprintf("%.2f", get_property(a, :occupancy, 1.0))
        bfac = @sprintf("%.2f", get_property(a, :tempfactor, 0.0))

        charge = if a.formal_charge == 0
            "?"
        elseif a.formal_charge > 0
            "$(a.formal_charge)+"
        else
            "$(-a.formal_charge)-"
        end

        model_num = a.frame_id

        println(io, "$group $( a.number) $type_sym $atom_name $alt_id $comp_id $chain_id $entity_id $seq_id $ins_code $x $y $z $occ $bfac $charge $seq_id $comp_id $chain_id $atom_name $model_num")
    end

    println(io, "#")
end

# ─── _struct_conn (disulfide bonds) ──────────────────────────────────

function _write_struct_conn(io::IO, ac::AbstractAtomContainer{T}) where T
    # Find disulfide bonds
    disulfides = filter(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(ac))
    isempty(disulfides) && return

    println(io, "loop_")
    println(io, "_struct_conn.id")
    println(io, "_struct_conn.conn_type_id")
    println(io, "_struct_conn.ptnr1_label_asym_id")
    println(io, "_struct_conn.ptnr1_label_comp_id")
    println(io, "_struct_conn.ptnr1_label_seq_id")
    println(io, "_struct_conn.ptnr1_label_atom_id")
    println(io, "_struct_conn.ptnr1_symmetry")
    println(io, "_struct_conn.ptnr2_label_asym_id")
    println(io, "_struct_conn.ptnr2_label_comp_id")
    println(io, "_struct_conn.ptnr2_label_seq_id")
    println(io, "_struct_conn.ptnr2_label_atom_id")
    println(io, "_struct_conn.ptnr2_symmetry")
    println(io, "_struct_conn.pdbx_dist_value")

    sys = parent_system(ac)

    for (i, b) in enumerate(disulfides)
        a1 = atom_by_idx(sys, b.a1)
        a2 = atom_by_idx(sys, b.a2)
        f1 = parent_fragment(a1)
        f2 = parent_fragment(a2)
        c1 = parent_chain(a1)
        c2 = parent_chain(a2)

        sym1 = get_property(b, :SYMMETRY_OPERATOR_0, 0)
        sym2 = get_property(b, :SYMMETRY_OPERATOR_1, 0)
        dist = get_property(b, :BOND_LENGTH, 0.0)

        sym1_str = sym1 == 0 ? "?" : string(sym1)
        sym2_str = sym2 == 0 ? "?" : string(sym2)
        dist_str = dist == 0.0 ? "?" : @sprintf("%.3f", dist)

        println(io, "disulf$(i) disulf $(c1.name) $(f1.name) $(f1.number) $(a1.name) $(sym1_str) $(c2.name) $(f2.name) $(f2.number) $(a2.name) $(sym2_str) $(dist_str)")
    end

    println(io, "#")
end

# ─── _struct_conf (helices) ──────────────────────────────────────────

function _write_struct_conf(io::IO, ac::AbstractAtomContainer{T}) where T
    helices = filter(ss -> ss.element == SecondaryStructureElement.Helix, secondary_structures(ac))
    isempty(helices) && return

    println(io, "loop_")
    println(io, "_struct_conf.conf_type_id")
    println(io, "_struct_conf.id")
    println(io, "_struct_conf.pdbx_PDB_helix_id")
    println(io, "_struct_conf.beg_label_comp_id")
    println(io, "_struct_conf.beg_label_asym_id")
    println(io, "_struct_conf.beg_label_seq_id")
    println(io, "_struct_conf.pdbx_beg_PDB_ins_code")
    println(io, "_struct_conf.end_label_comp_id")
    println(io, "_struct_conf.end_label_asym_id")
    println(io, "_struct_conf.end_label_seq_id")
    println(io, "_struct_conf.pdbx_end_PDB_ins_code")
    println(io, "_struct_conf.pdbx_PDB_helix_class")
    println(io, "_struct_conf.details")
    println(io, "_struct_conf.pdbx_PDB_helix_length")

    for ss in helices
        frags = fragments(ss)
        isempty(frags) && continue

        first_frag = first(frags)
        last_frag = last(frags)
        chain = parent_chain(ss)

        helix_class = get_property(ss, :HELIX_CLASS, 1)
        details = _cif_quote(get_property(ss, :COMMENT, "?"))
        helix_length = length(frags)

        beg_ins = get_property(first_frag, :insertion_code, "?")
        end_ins = get_property(last_frag, :insertion_code, "?")

        println(io, "HELX_P HELX_P$(ss.number) $(ss.name) $(first_frag.name) $(chain.name) $(first_frag.number) $(_cif_quote(beg_ins)) $(last_frag.name) $(chain.name) $(last_frag.number) $(_cif_quote(end_ins)) $(helix_class) $(details) $(helix_length)")
    end

    println(io, "#")
end

# ─── _struct_sheet_range (sheets) ────────────────────────────────────

function _write_struct_sheet_range(io::IO, ac::AbstractAtomContainer{T}) where T
    strands = filter(ss -> ss.element == SecondaryStructureElement.Strand, secondary_structures(ac))
    isempty(strands) && return

    println(io, "loop_")
    println(io, "_struct_sheet_range.sheet_id")
    println(io, "_struct_sheet_range.id")
    println(io, "_struct_sheet_range.beg_label_comp_id")
    println(io, "_struct_sheet_range.beg_label_asym_id")
    println(io, "_struct_sheet_range.beg_label_seq_id")
    println(io, "_struct_sheet_range.pdbx_beg_PDB_ins_code")
    println(io, "_struct_sheet_range.end_label_comp_id")
    println(io, "_struct_sheet_range.end_label_asym_id")
    println(io, "_struct_sheet_range.end_label_seq_id")
    println(io, "_struct_sheet_range.pdbx_end_PDB_ins_code")

    for ss in strands
        frags = fragments(ss)
        isempty(frags) && continue

        first_frag = first(frags)
        last_frag = last(frags)
        chain = parent_chain(ss)

        # Extract sheet name from "SheetName:RangeId" format used by PDB reader
        name_parts = split(ss.name, ":")
        sheet_id = length(name_parts) >= 1 ? name_parts[1] : "S1"

        beg_ins = get_property(first_frag, :insertion_code, "?")
        end_ins = get_property(last_frag, :insertion_code, "?")

        println(io, "$(sheet_id) $(ss.number) $(first_frag.name) $(chain.name) $(first_frag.number) $(_cif_quote(beg_ins)) $(last_frag.name) $(chain.name) $(last_frag.number) $(_cif_quote(end_ins))")
    end

    println(io, "#")
end
