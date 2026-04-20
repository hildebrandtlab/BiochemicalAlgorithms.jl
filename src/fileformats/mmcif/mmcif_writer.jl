using Printf

using BiochemicalAlgorithms

# Write an atom container as PDBx/mmCIF format to the given IO stream.
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

# Quote a string value for inline CIF output (a single token within a row).
# Always returns a single-line representation: a CIF semicolon text block
# must start at column 1 on its own line, and our caller `join`s fields with
# spaces, so we can't emit one here. Newlines and carriage returns are
# collapsed to spaces; embedded quotes are handled by switching delimiters.
function _cif_quote(s::AbstractString)
    isempty(s) && return "."
    if occursin('\n', s) || occursin('\r', s)
        s = replace(s, r"[\r\n]+" => " ")
    end
    # No quoting needed for simple values.
    if !any(c -> isspace(c), s) && !startswith(s, '_') && !startswith(s, '#') &&
       !startswith(s, '\'') && !startswith(s, '"') && s != "." && s != "?"
        return s
    end
    # Use single quotes if possible.
    if !occursin('\'', s)
        return "'$s'"
    end
    # Use double quotes.
    if !occursin('"', s)
        return "\"$s\""
    end
    # Both quote types present — escape internal single quotes and wrap in
    # single quotes. Not strictly CIF-compliant (CIF has no escape syntax),
    # but produces parseable output for the round-trip case where neither
    # quote alone fits.
    return "'" * replace(s, '\'' => "\\'") * "'"
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
        # mmCIF requires `_atom_site.type_symbol` to be a chemical symbol or
        # a missing-value marker. Map `Elements.Unknown` to `?` rather than
        # the literal string `"Unknown"`. Preserve deuterium as `D` when the
        # `:is_deuterium` flag is set (the reader maps `D` to `Elements.H`).
        type_sym = if has_flag(a, :is_deuterium)
            "D"
        elseif a.element == Elements.Unknown
            "?"
        else
            string(a.element)
        end
        atom_name = _cif_quote(a.name)
        alt_id = "."

        comp_id  = _cif_quote(isnothing(frag)  ? "UNK" : frag.name)
        chain_id = _cif_quote(isnothing(chain) ? "."   : chain.name)
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

        println(io, join((
            group, a.number, type_sym, atom_name, alt_id,
            comp_id, chain_id, entity_id, seq_id, ins_code,
            x, y, z, occ, bfac, charge,
            seq_id, comp_id, chain_id, atom_name, model_num,
        ), " "))
    end

    println(io, "#")
end

# ─── _struct_conn (inter-residue bonds) ──────────────────────────────

# Map BondOrder to mmCIF `pdbx_value_order` strings.
const _BOND_ORDER_TO_CIF = Dict(
    BondOrder.Single    => "sing",
    BondOrder.Double    => "doub",
    BondOrder.Triple    => "trip",
    BondOrder.Quadruple => "quad",
    BondOrder.Aromatic  => "arom",
)

# Pick the mmCIF conn_type_id for a bond. Prefer the original parsed value
# (stored as :CIF_CONN_TYPE on read) so reader→writer round-trips don't lose
# the covale/metalc distinction. Fall back to the bond's flags otherwise.
function _conn_type_for(b)
    t = get_property(b, :CIF_CONN_TYPE, nothing)
    !isnothing(t) && return String(t)
    has_flag(b, :TYPE__DISULPHIDE_BOND) && return "disulf"
    has_flag(b, :TYPE__HYDROGEN)        && return "hydrog"
    has_flag(b, :TYPE__SALT_BRIDGE)     && return "saltbr"
    has_flag(b, :TYPE__COVALENT)        && return "covale"
    return nothing
end

function _write_struct_conn(io::IO, ac::AbstractAtomContainer{T}) where T
    # Bonds we emit: anything that has a known mmCIF conn_type. Programmatically
    # constructed bonds without flags are skipped (they're not inter-residue
    # connections in the mmCIF sense).
    conn_bonds = [(b, _conn_type_for(b)) for b in bonds(ac)]
    filter!(p -> !isnothing(p[2]), conn_bonds)
    isempty(conn_bonds) && return

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
    println(io, "_struct_conn.pdbx_value_order")

    sys = parent_system(ac)

    for (i, (b, ctype)) in enumerate(conn_bonds)
        a1, a2 = get_partners(b)
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
        order_str = get(_BOND_ORDER_TO_CIF, b.order, "?")

        println(io, join((
            "$(ctype)$i", ctype,
            _cif_quote(c1.name), _cif_quote(f1.name), f1.number, _cif_quote(a1.name), sym1_str,
            _cif_quote(c2.name), _cif_quote(f2.name), f2.number, _cif_quote(a2.name), sym2_str,
            dist_str, order_str,
        ), " "))
    end

    println(io, "#")
end

# ─── _struct_conf (helices) ──────────────────────────────────────────

function _write_struct_conf(io::IO, ac::AbstractAtomContainer{T}) where T
    helices = filter(ss -> ss.type == SecondaryStructureElement.Helix, secondary_structures(ac))
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

        beg_ins = _cif_quote(get_property(first_frag, :insertion_code, "?"))
        end_ins = _cif_quote(get_property(last_frag, :insertion_code, "?"))

        println(io, join((
            "HELX_P", "HELX_P$(ss.number)", _cif_quote(ss.name),
            _cif_quote(first_frag.name), _cif_quote(chain.name), first_frag.number, beg_ins,
            _cif_quote(last_frag.name),  _cif_quote(chain.name), last_frag.number,  end_ins,
            helix_class, details, helix_length,
        ), " "))
    end

    println(io, "#")
end

# ─── _struct_sheet_range (sheets) ────────────────────────────────────

function _write_struct_sheet_range(io::IO, ac::AbstractAtomContainer{T}) where T
    strands = filter(ss -> ss.type == SecondaryStructureElement.Strand, secondary_structures(ac))
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

        # Recover the original sheet_id and range_id from the name. The PDB
        # reader's postprocessing stores the strand name as `<sheet>:<range>`
        # but renumbers `ss.number` sequentially within the chain, so the
        # original range_id is only available via the name.
        name_parts = split(ss.name, ":")
        sheet_id = first(name_parts)
        range_id = length(name_parts) >= 2 ? something(tryparse(Int, name_parts[2]), ss.number) : ss.number

        beg_ins = _cif_quote(get_property(first_frag, :insertion_code, "?"))
        end_ins = _cif_quote(get_property(last_frag, :insertion_code, "?"))

        println(io, join((
            _cif_quote(sheet_id), range_id,
            _cif_quote(first_frag.name), _cif_quote(chain.name), first_frag.number, beg_ins,
            _cif_quote(last_frag.name),  _cif_quote(chain.name), last_frag.number,  end_ins,
        ), " "))
    end

    println(io, "#")
end
