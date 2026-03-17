using BiochemicalAlgorithms
using BiochemicalAlgorithms: CIFFile, CIFDataBlock, CIFLoop, CIFParser,
    parse_cif_file, parse_cif, parse_element_string, _fragment_variant
import BiochemicalAlgorithms: PDBDetails
using DataStructures: Deque

"""
    read_mmcif(fname_io::Union{AbstractString, IO}, ::Type{T} = Float32; create_coils::Bool = true) -> System{T}

Read a PDBx/mmCIF file and return a System.

Models are stored as frames, using the model number as `frame_id`.
"""
function read_mmcif(fname_io::Union{AbstractString, IO}, ::Type{T} = Float32;
        create_coils::Bool = true) where {T <: Real}
    cif = if fname_io isa AbstractString
        parse_cif_file(fname_io)
    else
        parse_cif(fname_io)
    end

    if isempty(cif.blocks)
        error("mmCIF file contains no data blocks")
    end

    block = first(values(cif.blocks))

    # Set system name from data block name (or filename if available)
    sys_name = if fname_io isa AbstractString
        bn = basename(fname_io)
        # strip extension
        idx = findlast('.', bn)
        isnothing(idx) ? bn : bn[1:idx-1]
    else
        block.name
    end

    sys = System{T}(sys_name)
    mol = Molecule(sys; name = sys_name)

    # Parse _atom_site loop
    atom_loop = _find_loop(block, "_atom_site.")
    if isnothing(atom_loop)
        return sys
    end

    _build_atoms!(sys, mol, atom_loop, T)

    # Build fragment cache for postprocessing
    fragment_cache = Dict{PDBDetails.UniqueResidueID, Fragment{T}}()
    for f in fragments(sys)
        fragment_cache[PDBDetails.UniqueResidueID(
            f.name,
            parent_chain(f).name,
            f.number,
            get_property(f, :insertion_code, " ")
        )] = f
    end

    # Parse disulfide bonds from _struct_conn
    ssbonds = _parse_ssbonds(block)
    PDBDetails.postprocess_ssbonds_!(sys, ssbonds, fragment_cache)

    # Parse secondary structure from _struct_conf and _struct_sheet_range
    ss_records = _parse_secondary_structures(block)
    PDBDetails.postprocess_secondary_structures_!(sys, ss_records, fragment_cache, create_coils)

    sys
end

# ─── Helpers ──────────────────────────────────────────────────────────

"""Find a loop in the data block whose tags start with the given prefix."""
function _find_loop(block::CIFDataBlock, prefix::String)
    for loop in block.loops
        if !isempty(loop.tags) && startswith(loop.tags[1], prefix)
            return loop
        end
    end
    return nothing
end

"""Build a tag→column-index map for a loop."""
@inline function _col_map(loop::CIFLoop)
    Dict(tag => i for (i, tag) in enumerate(loop.tags))
end

"""Get a value from a row, returning `nothing` for CIF missing values (? and .)."""
@inline function _get(row::Vector{String}, col::Int)
    v = row[col]
    (v == "?" || v == ".") ? nothing : v
end

"""Get a value with a default for missing."""
@inline function _get(row::Vector{String}, col::Int, default)
    v = _get(row, col)
    isnothing(v) ? default : v
end

"""Get an optional column index, returning nothing if the tag is not present."""
@inline function _optcol(cols::Dict{String, Int}, tag::String)
    get(cols, tag, nothing)
end

# ─── Atom Site Parsing ────────────────────────────────────────────────

function _build_atoms!(sys::System{T}, mol::Molecule{T}, loop::CIFLoop, ::Type{T}) where {T <: Real}
    cols = _col_map(loop)

    # Required columns
    c_group     = cols["_atom_site.group_PDB"]
    c_id        = cols["_atom_site.id"]
    c_symbol    = cols["_atom_site.type_symbol"]
    c_atom_name = cols["_atom_site.label_atom_id"]
    c_alt_id    = cols["_atom_site.label_alt_id"]
    c_comp_id   = cols["_atom_site.label_comp_id"]
    c_asym_id   = cols["_atom_site.label_asym_id"]
    c_seq_id    = cols["_atom_site.label_seq_id"]
    c_x         = cols["_atom_site.Cartn_x"]
    c_y         = cols["_atom_site.Cartn_y"]
    c_z         = cols["_atom_site.Cartn_z"]

    # Optional columns (use auth_* when available, fall back to label_*)
    c_auth_asym  = _optcol(cols, "_atom_site.auth_asym_id")
    c_auth_comp  = _optcol(cols, "_atom_site.auth_comp_id")
    c_auth_seq   = _optcol(cols, "_atom_site.auth_seq_id")
    c_auth_atom  = _optcol(cols, "_atom_site.auth_atom_id")
    c_ins_code   = _optcol(cols, "_atom_site.pdbx_PDB_ins_code")
    c_occupancy  = _optcol(cols, "_atom_site.occupancy")
    c_bfactor    = _optcol(cols, "_atom_site.B_iso_or_equiv")
    c_charge     = _optcol(cols, "_atom_site.pdbx_formal_charge")
    c_model      = _optcol(cols, "_atom_site.pdbx_PDB_model_num")

    current_chain::Union{Chain{T}, Nothing} = nothing
    current_frag::Union{Fragment{T}, Nothing} = nothing
    current_chain_id = ""
    current_frag_key = ("", 0, "")  # (comp_id, seq_id, ins_code)

    altloc_warning = false
    first_altloc_id = Dict{Tuple{String, Int, String, String}, String}()  # (chain, seq, ins, atom_name) -> first altloc

    for row in loop.rows
        # Alternate location handling
        alt_id_raw = row[c_alt_id]
        alt_id = (alt_id_raw == "." || alt_id_raw == "?") ? nothing : alt_id_raw

        # Use auth_* fields when available (matches PDB convention)
        chain_id = isnothing(c_auth_asym) ? row[c_asym_id] : _get(row, c_auth_asym, row[c_asym_id])
        comp_id  = isnothing(c_auth_comp) ? row[c_comp_id] : _get(row, c_auth_comp, row[c_comp_id])
        atom_name = isnothing(c_auth_atom) ? row[c_atom_name] : _get(row, c_auth_atom, row[c_atom_name])

        seq_id_str = isnothing(c_auth_seq) ? _get(row, c_seq_id, "0") : _get(row, c_auth_seq, _get(row, c_seq_id, "0"))
        seq_id = tryparse(Int, seq_id_str)
        if isnothing(seq_id)
            seq_id = 0
        end

        ins_code = isnothing(c_ins_code) ? " " : _get(row, c_ins_code, " ")

        # Skip alternate locations beyond the first
        if !isnothing(alt_id)
            key = (chain_id, seq_id, ins_code, atom_name)
            existing = get(first_altloc_id, key, nothing)
            if isnothing(existing)
                first_altloc_id[key] = alt_id
            elseif alt_id != existing
                if !altloc_warning
                    @warn "load_mmcif: alternate locations other than $(existing) are currently not supported. Affected records have been ignored!"
                    altloc_warning = true
                end
                continue
            end
        end

        # Chain
        if isnothing(current_chain) || chain_id != current_chain_id
            current_chain = Chain(mol; name = chain_id)
            current_chain_id = chain_id
            current_frag = nothing
            current_frag_key = ("", 0, "")
        end

        # Fragment
        frag_key = (comp_id, seq_id, ins_code)
        if isnothing(current_frag) || frag_key != current_frag_key
            is_hetero = strip(row[c_group]) == "HETATM"
            current_frag = Fragment(current_chain, seq_id;
                name = comp_id,
                variant = _fragment_variant(comp_id),
                properties = Properties([
                    :is_hetero_fragment => is_hetero,
                    :insertion_code => ins_code
                ])
            )
            current_frag_key = frag_key
        end

        # Atom
        serial = parse(Int, row[c_id])
        element = parse_element_string(strip(row[c_symbol]))

        x = parse(T, row[c_x])
        y = parse(T, row[c_y])
        z = parse(T, row[c_z])

        occupancy = isnothing(c_occupancy) ? T(1.0) : T(parse(Float64, _get(row, c_occupancy, "1.0")))
        tempfactor = isnothing(c_bfactor) ? T(0.0) : T(parse(Float64, _get(row, c_bfactor, "0.0")))

        charge_str = isnothing(c_charge) ? nothing : _get(row, c_charge)
        formal_charge = isnothing(charge_str) ? Int(0) : (tryparse(Int, charge_str) === nothing ? 0 : parse(Int, charge_str))

        model_num = isnothing(c_model) ? 1 : parse(Int, _get(row, c_model, "1"))

        is_hetero_atom = strip(row[c_group]) == "HETATM"
        is_deuterium = strip(row[c_symbol]) == "D"

        flags = Flags()
        is_hetero_atom && push!(flags, :is_hetero_atom)
        is_deuterium && push!(flags, :is_deuterium)

        Atom(current_frag, serial, element;
            name = atom_name,
            r = Vector3{T}(x, y, z),
            formal_charge = formal_charge,
            properties = Properties([
                :tempfactor => tempfactor,
                :occupancy => occupancy,
                :is_hetero_atom => is_hetero_atom,
                :insertion_code => ins_code
            ]),
            flags = flags,
            frame_id = model_num
        )
    end
end

# ─── SSBond Parsing ──────────────────────────────────────────────────

function _parse_ssbonds(block::CIFDataBlock)
    ssbonds = Deque{PDBDetails.SSBondRecord}()
    loop = _find_loop(block, "_struct_conn.")
    isnothing(loop) && return ssbonds

    cols = _col_map(loop)

    c_type = cols["_struct_conn.conn_type_id"]

    # Use auth fields for residue identification
    c_p1_asym = get(cols, "_struct_conn.ptnr1_auth_asym_id", get(cols, "_struct_conn.ptnr1_label_asym_id", 0))
    c_p1_comp = get(cols, "_struct_conn.ptnr1_auth_comp_id", get(cols, "_struct_conn.ptnr1_label_comp_id", 0))
    c_p1_seq  = get(cols, "_struct_conn.ptnr1_auth_seq_id", get(cols, "_struct_conn.ptnr1_label_seq_id", 0))
    c_p2_asym = get(cols, "_struct_conn.ptnr2_auth_asym_id", get(cols, "_struct_conn.ptnr2_label_asym_id", 0))
    c_p2_comp = get(cols, "_struct_conn.ptnr2_auth_comp_id", get(cols, "_struct_conn.ptnr2_label_comp_id", 0))
    c_p2_seq  = get(cols, "_struct_conn.ptnr2_auth_seq_id", get(cols, "_struct_conn.ptnr2_label_seq_id", 0))

    c_p1_ins  = _optcol(cols, "_struct_conn.pdbx_ptnr1_PDB_ins_code")
    c_p2_ins  = _optcol(cols, "_struct_conn.pdbx_ptnr2_PDB_ins_code")
    c_sym1    = _optcol(cols, "_struct_conn.ptnr1_symmetry")
    c_sym2    = _optcol(cols, "_struct_conn.ptnr2_symmetry")
    c_dist    = _optcol(cols, "_struct_conn.pdbx_dist_value")

    n = 0
    for row in loop.rows
        strip(row[c_type]) == "disulf" || continue
        n += 1

        p1_ins = isnothing(c_p1_ins) ? " " : _get(row, c_p1_ins, " ")
        p2_ins = isnothing(c_p2_ins) ? " " : _get(row, c_p2_ins, " ")

        first_res = PDBDetails.UniqueResidueID(
            strip(row[c_p1_comp]),
            strip(row[c_p1_asym]),
            parse(Int, strip(row[c_p1_seq])),
            p1_ins
        )
        second_res = PDBDetails.UniqueResidueID(
            strip(row[c_p2_comp]),
            strip(row[c_p2_asym]),
            parse(Int, strip(row[c_p2_seq])),
            p2_ins
        )

        sym1 = isnothing(c_sym1) ? 0 : _parse_symmetry_operator(_get(row, c_sym1, "0"))
        sym2 = isnothing(c_sym2) ? 0 : _parse_symmetry_operator(_get(row, c_sym2, "0"))
        dist = isnothing(c_dist) ? 0.0 : parse(Float64, _get(row, c_dist, "0.0"))

        push!(ssbonds, PDBDetails.SSBondRecord(n, first_res, second_res, sym1, sym2, dist))
    end

    return ssbonds
end

"""Parse CIF symmetry operator string like '1_555' into an integer."""
function _parse_symmetry_operator(s::String)
    s = strip(s)
    (s == "?" || s == "." || isempty(s)) && return 0
    # Try parsing as integer directly
    v = tryparse(Int, replace(s, "_" => ""))
    isnothing(v) ? 0 : v
end

# ─── Secondary Structure Parsing ─────────────────────────────────────

function _parse_secondary_structures(block::CIFDataBlock)
    records = Deque{Union{PDBDetails.HelixRecord, PDBDetails.SheetRecord, PDBDetails.TurnRecord}}()

    _parse_helices!(records, block)
    _parse_sheets!(records, block)

    return records
end

function _parse_helices!(records, block::CIFDataBlock)
    loop = _find_loop(block, "_struct_conf.")
    isnothing(loop) && return

    cols = _col_map(loop)

    c_id = cols["_struct_conf.id"]

    # Use auth fields when available
    c_beg_comp = get(cols, "_struct_conf.beg_auth_comp_id", get(cols, "_struct_conf.beg_label_comp_id", 0))
    c_beg_asym = get(cols, "_struct_conf.beg_auth_asym_id", get(cols, "_struct_conf.beg_label_asym_id", 0))
    c_beg_seq  = get(cols, "_struct_conf.beg_auth_seq_id", get(cols, "_struct_conf.beg_label_seq_id", 0))
    c_end_comp = get(cols, "_struct_conf.end_auth_comp_id", get(cols, "_struct_conf.end_label_comp_id", 0))
    c_end_asym = get(cols, "_struct_conf.end_auth_asym_id", get(cols, "_struct_conf.end_label_asym_id", 0))
    c_end_seq  = get(cols, "_struct_conf.end_auth_seq_id", get(cols, "_struct_conf.end_label_seq_id", 0))

    c_helix_id    = _optcol(cols, "_struct_conf.pdbx_PDB_helix_id")
    c_helix_class = _optcol(cols, "_struct_conf.pdbx_PDB_helix_class")
    c_details     = _optcol(cols, "_struct_conf.details")
    c_beg_ins     = _optcol(cols, "_struct_conf.pdbx_beg_PDB_ins_code")
    c_end_ins     = _optcol(cols, "_struct_conf.pdbx_end_PDB_ins_code")

    n = 0
    for row in loop.rows
        n += 1

        helix_name = isnothing(c_helix_id) ? _get(row, c_id, "H$n") : _get(row, c_helix_id, "H$n")
        beg_ins = isnothing(c_beg_ins) ? " " : _get(row, c_beg_ins, " ")
        end_ins = isnothing(c_end_ins) ? " " : _get(row, c_end_ins, " ")

        initial = PDBDetails.UniqueResidueID(
            strip(row[c_beg_comp]),
            strip(row[c_beg_asym]),
            parse(Int, strip(row[c_beg_seq])),
            beg_ins
        )
        terminal = PDBDetails.UniqueResidueID(
            strip(row[c_end_comp]),
            strip(row[c_end_asym]),
            parse(Int, strip(row[c_end_seq])),
            end_ins
        )

        helix_class = isnothing(c_helix_class) ? 1 : (tryparse(Int, _get(row, c_helix_class, "1")) === nothing ? 1 : parse(Int, _get(row, c_helix_class, "1")))
        details = isnothing(c_details) ? "" : _get(row, c_details, "")

        push!(records, PDBDetails.HelixRecord(n, helix_name, initial, terminal, helix_class, details))
    end
end

function _parse_sheets!(records, block::CIFDataBlock)
    loop = _find_loop(block, "_struct_sheet_range.")
    isnothing(loop) && return

    cols = _col_map(loop)

    c_sheet_id = cols["_struct_sheet_range.sheet_id"]
    c_range_id = cols["_struct_sheet_range.id"]

    # Use auth fields when available
    c_beg_comp = get(cols, "_struct_sheet_range.beg_auth_comp_id", get(cols, "_struct_sheet_range.beg_label_comp_id", 0))
    c_beg_asym = get(cols, "_struct_sheet_range.beg_auth_asym_id", get(cols, "_struct_sheet_range.beg_label_asym_id", 0))
    c_beg_seq  = get(cols, "_struct_sheet_range.beg_auth_seq_id", get(cols, "_struct_sheet_range.beg_label_seq_id", 0))
    c_end_comp = get(cols, "_struct_sheet_range.end_auth_comp_id", get(cols, "_struct_sheet_range.end_label_comp_id", 0))
    c_end_asym = get(cols, "_struct_sheet_range.end_auth_asym_id", get(cols, "_struct_sheet_range.end_label_asym_id", 0))
    c_end_seq  = get(cols, "_struct_sheet_range.end_auth_seq_id", get(cols, "_struct_sheet_range.end_label_seq_id", 0))

    c_beg_ins = _optcol(cols, "_struct_sheet_range.pdbx_beg_PDB_ins_code")
    c_end_ins = _optcol(cols, "_struct_sheet_range.pdbx_end_PDB_ins_code")

    # Try to get sense from _struct_sheet_order
    sense_map = _parse_sheet_order_sense(block)

    for row in loop.rows
        sheet_id = strip(row[c_sheet_id])
        range_id = parse(Int, strip(row[c_range_id]))

        beg_ins = isnothing(c_beg_ins) ? " " : _get(row, c_beg_ins, " ")
        end_ins = isnothing(c_end_ins) ? " " : _get(row, c_end_ins, " ")

        initial = PDBDetails.UniqueResidueID(
            strip(row[c_beg_comp]),
            strip(row[c_beg_asym]),
            parse(Int, strip(row[c_beg_seq])),
            beg_ins
        )
        terminal = PDBDetails.UniqueResidueID(
            strip(row[c_end_comp]),
            strip(row[c_end_asym]),
            parse(Int, strip(row[c_end_seq])),
            end_ins
        )

        sense = get(sense_map, (sheet_id, range_id), 0)

        push!(records, PDBDetails.SheetRecord(range_id, sheet_id, initial, terminal, sense))
    end
end

"""Parse _struct_sheet_order to get sense values for each sheet range."""
function _parse_sheet_order_sense(block::CIFDataBlock)
    sense_map = Dict{Tuple{String, Int}, Int}()
    loop = _find_loop(block, "_struct_sheet_order.")
    isnothing(loop) && return sense_map

    cols = _col_map(loop)
    c_sheet = cols["_struct_sheet_order.sheet_id"]
    c_range2 = cols["_struct_sheet_order.range_id_2"]
    c_sense = _optcol(cols, "_struct_sheet_order.sense")

    isnothing(c_sense) && return sense_map

    for row in loop.rows
        sheet_id = strip(row[c_sheet])
        range_id = parse(Int, strip(row[c_range2]))
        sense_str = _get(row, c_sense, "parallel")
        sense = sense_str == "anti-parallel" ? -1 : (sense_str == "parallel" ? 1 : 0)
        sense_map[(sheet_id, range_id)] = sense
    end

    return sense_map
end
