using Printf
using Dates
using DataFrames

using BiochemicalAlgorithms

function extract_structure(ac::AtomContainer)
    # TODO: handle crystal information, like in the original BALL PDB writer

    # Walk over all atoms and collect all chains, bonds, etc...


end

# this really should be in the standard library
unval(::Val{x}) where x = x

function fix_nucleotide_name(name::String)
    replace(name, 
            r"^DA$" => "A",
            r"^DC$" => "C",
            r"^DG$" => "G",
            r"^DT$" => "T",
            r"^DI$" => "I",
            r"^DU$" => "U")
end

function residue_details(res::Fragment{T}) where {T <:Real}
    res.name, isnothing(parent_chain(res)) ? "" : parent_chain(res).name, res.number, get_property(res, :insertion_code, " ")
end

function atom_details(atom::Atom{T}) where {T <:Real}
    res = parent_fragment(atom)

    if isnothing(res)
        "UNK", "", 0, ""
    else
        residue_details(res)
    end
end

function extract_records(record_type::Val, pdb_info::PDBInfo)
    records = Iterators.filter(r -> Val(r.type) == record_type, pdb_info.records)
end

function break_into_continuation_lines(line::String, c::Int = 1)
    if length(line) > 70 # we need to leave room for the eol, so we skip one character at the end
        [(c, line[1:69]), break_into_continuation_lines(String(line[71:end]), c+1)...]
    else
        [(c, line)]
    end
end

function write_record(io::IO, pdb_info::PDBInfo, tag_name::String, args...)
    type = RECORD_MAP[tag_name]
    
    # CONECT-records are handled differently
    if type.record_type === RECORD_TYPE__REMARK
        pdb_info.writer_stats.remark_records += 1
    elseif type.record_type === RECORD_TYPE__HET
        pdb_info.writer_stats.het_records += 1
    elseif type.record_type === RECORD_TYPE__HELIX
        pdb_info.writer_stats.helix_records += 1
    elseif type.record_type === RECORD_TYPE__SHEET
        pdb_info.writer_stats.sheet_records += 1
    elseif type.record_type === RECORD_TYPE__TURN
        pdb_info.writer_stats.turn_records += 1
    elseif type.record_type === RECORD_TYPE__SITE
        pdb_info.writer_stats.site_records += 1
    elseif type.record_type in [
            RECORD_TYPE__ORIGX1, RECORD_TYPE__ORIGX2, RECORD_TYPE__ORIGX3,
            RECORD_TYPE__SCALE1, RECORD_TYPE__SCALE2, RECORD_TYPE__SCALE3,
            RECORD_TYPE__MTRIX1, RECORD_TYPE__MTRIX2, RECORD_TYPE__MTRIX3    
        ]
        pdb_info.writer_stats.coordinate_transformation_records += 1
    elseif type.record_type in [RECORD_TYPE__ATOM, RECORD_TYPE__HETATM]
        pdb_info.writer_stats.atomic_coordinate_records += 1
    elseif type.record_type == RECORD_TYPE__TER
        pdb_info.writer_stats.ter_records += 1
    else type.record_type == RECORD_TYPE__SEQRES
        pdb_info.writer_stats.seqres_records += 1
    end

    Printf.format(io, Printf.Format(type.tag * type.format_string * "\n"), args...)
end

function write_record(io::IO, pdb_info::PDBInfo, ::Val{RECORD_TYPE__TITLE}, continuation::Int, title_line::String)
    if continuation < 2
        Printf.format(io, Printf.Format(RECORD_TAG_TITLE * "    %-70.70s\n"), title_line)
    else
        write_record(io, pdb_info, RECORD_TAG_TITLE, string(continuation), title_line)
    end
end

function write_record(io::IO, pdb_info::PDBInfo, ::Val{RECORD_TYPE__END}, data::Tuple)
    write_record(io, pdb_info, RECORD_TAG_END, data...)
end

function write_record(io::IO, pdb_info::PDBInfo, ::Val{RECORD_TYPE__SHEET}, data...)
    pdb_info.writer_stats.sheet_records += 1

    # we don't store all the data required for the SHEET record, so
    # we have to massage the format string here
    fstring = Printf.Format(RECORD_TAG_SHEET * join(split(FORMAT_SHEET, "%")[1:end-10], "%") * repeat(" ", 39) *"\n")
    Printf.format(io, fstring, data...)
end

function write_record(io::IO, pdb_info::PDBInfo, ::Val{RECORD_TYPE__CONECT}, data...)
    pdb_info.writer_stats.conect_records += 1

    # we generally use the 2006 format for the CONECT record, but we have to massage the format string
    # because we might have too few indices if there are less than 4 bonds for this atom
    fstring = Printf.Format(RECORD_TAG_CONECT * replace(FORMAT_CON064, r"%(\d+)ld" => s"%\1s") * "\n")
    Printf.format(io, fstring, data...)
end

function write_record(io::IO, pdb_info::PDBInfo, record_type::Val, data...)
    write_record(io, pdb_info, typeformat_by_type(unval(record_type)).tag, data...)
end

function write_record(io::IO, pdb_info::PDBInfo, record::PDBRecord)
    write_record(io, pdb_info, Val(record.type), record.data...)
end

function write_records(io::IO, pdb_info::PDBInfo, ::Val{RECORD_TYPE__TITLE})
    title_lines = split(pdb_info.title, "\n")

    for line in title_lines
        for (c, line) in break_into_continuation_lines(String(line))
            write_record(io, pdb_info, Val(RECORD_TYPE__TITLE), c, line)
        end
    end
end

function write_records(io::IO, pdb_info::PDBInfo, record_type)
    records = extract_records(record_type, pdb_info)

    for record in records
        write_record(io, pdb_info, record)
    end
end

function write_seqres(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # iterate over all chains
    for chain in chains(ac)
        res = fragments(chain)
        nres = length(res)
        # each chain is stored in groups of 13 residues
        for (i, rs) in enumerate(Iterators.partition(res, 13))
            rs = vcat(map(r -> fix_nucleotide_name(r.name), rs), repeat([""], 13 - length(rs)))
            write_record(io, pdb_info, RECORD_TAG_SEQRES, i, chain.name, nres, rs...)
        end
    end
end

function write_title_section(io::IO, pdb_info::PDBInfo)
    # --- HEADER ---

    # Retrieve the current date...
    date = pdb_info.deposition_date != "" ? pdb_info.deposition_date : uppercase(Dates.format(now(), "dd-u-yy"))
    write_record(io, pdb_info, "HEADER", pdb_info.name, date, pdb_info.id)

    # --- OBSLTE ---
    write_records(io, pdb_info, Val(RECORD_TYPE__OBSLTE))
    # --- TITLE ---
    write_records(io, pdb_info, Val(RECORD_TYPE__TITLE))
    # --- CAVEAT ---
    write_records(io, pdb_info, Val(RECORD_TYPE__CAVEAT))
    # --- COMPND ---
    write_records(io, pdb_info, Val(RECORD_TYPE__COMPND))
    # --- SOURCE ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SOURCE))
    # --- KEYWDS ---
    write_records(io, pdb_info, Val(RECORD_TYPE__KEYWDS))
    # --- EXPDTA ---
    write_records(io, pdb_info, Val(RECORD_TYPE__EXPDTA))
    # --- AUTHOR ---
    write_records(io, pdb_info, Val(RECORD_TYPE__AUTHOR))
    # --- REVDAT ---
    write_records(io, pdb_info, Val(RECORD_TYPE__REVDAT))
    # --- SPRSDE ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SPRSDE))
    # --- JRNL ---
    write_records(io, pdb_info, Val(RECORD_TYPE__JRNL))
    # --- REMARK ---
    write_records(io, pdb_info, Val(RECORD_TYPE__REMARK))
end

function write_primary_structure_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # --- DBREF ---
    write_records(io, pdb_info, Val(RECORD_TYPE__DBREF))
    # --- SEQADV ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SEQADV))
    # --- SEQRES ---
    write_seqres(io, pdb_info, ac)
    # --- MODRES ---
    write_records(io, pdb_info, Val(RECORD_TYPE__MODRES))
end

function write_heterogen_section(io::IO, pdb_info::PDBInfo)
    # --- HET ---
    write_records(io, pdb_info, Val(RECORD_TYPE__HET))
    # --- HETNAM ---
    write_records(io, pdb_info, Val(RECORD_TYPE__HETNAM))
    # --- HETSYN ---
    write_records(io, pdb_info, Val(RECORD_TYPE__HETSYN))
    # --- FORMUL ---
    write_records(io, pdb_info, Val(RECORD_TYPE__FORMUL))
end

function write_helix_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    helices = filter(
        ss -> ss.type == SecondaryStructureElement.Helix, 
        secondary_structures(ac)
    )

    for (helix_number, helix) in enumerate(helices)
        rs = fragments(helix)

        helix_start = first(rs)
        helix_end   =  last(rs)

        helix_class   = get_property(helix, :HELIX_CLASS, 1)
        helix_comment = get_property(helix, :COMMENT,    "")

        write_record(io, pdb_info, RECORD_TAG_HELIX, 
            helix_number, helix.name, 
            residue_details(helix_start)...,
            residue_details(helix_end)...,
            helix_class, helix_comment,
            length(rs)
        )
    end
end

function write_sheet_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # this uses the same convention as BALL did, but we might have to 
    # check this eventually...
    # this is the original BALL comment:
    #   First, count the number of strands. Not exactly what 
	#   we need, but it'll do for now. We would have to count
	#   the strands per sheet for every distinct sheet ID and
	#   account for barrels (closed sheets: last strand = first strand)
	#   as well.

    strands = filter(
        ss -> ss.type == SecondaryStructureElement.Strand, 
        secondary_structures(ac)
    )

    for strand in strands
        rs = fragments(strand)

        name_parts = split(strand.name, ":")

        strand_number = 
            if length(name_parts) == 1
                1
            else
                try
                    parse(Int, name_parts[2])
                catch
                    1
                end
            end
        strand_name = name_parts[1]

        strand_start = first(rs)
        strand_end   =  last(rs)

        strand_sense = get_property(strand, :STRAND_SENSE, 0)

        write_record(io, pdb_info, Val(RECORD_TYPE__SHEET),
            strand_number, strand_name, length(strands),
            residue_details(strand_start)...,
            residue_details(strand_end)...,
            strand_sense
        )
    end
end

function write_turn_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    turns = filter(
        ss -> ss.type == SecondaryStructureElement.Turn, 
        secondary_structures(ac)
    )

    for (turn_number, turn) in enumerate(turns)
        write(io, pdb_info, RECORD_TAG_TURN, 
            turn_number, turn.name,
            residue_details(first(turn))...,
            residue_details(last(turn))...,
            get_property(turn, :COMMENT, "")
        )
    end
end

function write_secondary_structure_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # --- HELIX ---
    write_helix_section(io, pdb_info, ac)
    # --- SHEET ---
    write_sheet_section(io, pdb_info, ac)
    # --- TURN ---
    write_turn_section(io, pdb_info, ac)
end

function write_ssbond_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    ssbonds = filter(
        ss -> has_flag(ss, :TYPE__DISULPHIDE_BOND),
        bonds(ac)
    )

    for (ssbond_number, ssbond) in enumerate(ssbonds)
        first_partner, second_partner = get_partners(ssbond)
        symmetry_operator_0 = get_property(ssbond, :SYMMETRY_OPERATOR_0, 0)
        symmetry_operator_1 = get_property(ssbond, :SYMMETRY_OPERATOR_1, 0)
        bond_length = get_property(ssbond, :BOND_LENGTH, 0.0)

        write_record(io, pdb_info, RECORD_TAG_SSBOND, 
            ssbond_number, 
            atom_details(first_partner)...,
            atom_details(second_partner)...,
            symmetry_operator_0, symmetry_operator_1, bond_length
        )
    end
end

function write_connectivity_annotation_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # --- SSBOND ---
    write_ssbond_section(io, pdb_info, ac)
    # --- LINK ---
    write_records(io, pdb_info, Val(RECORD_TYPE__LINK))
    # --- HYDBND ---
    write_records(io, pdb_info, Val(RECORD_TYPE__HYDBND))
    # --- SLTBRG ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SLTBRG))
    # --- CISPEP ---
    write_records(io, pdb_info, Val(RECORD_TYPE__CISPEP))
end

function write_miscellaneous_features_section(io::IO, pdb_info::PDBInfo)
    # --- SITE ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SITE))
end

function write_crystallographic_section(io::IO, pdb_info::PDBInfo)
    # --- CRYST1 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__CRYST1))
    # --- ORIGX1 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__ORIGX1))
    # --- ORIGX2 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__ORIGX2))
    # --- ORIGX3 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__ORIGX3))
    # --- SCALE1 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SCALE1))
    # --- SCALE2 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SCALE2))
    # --- SCALE3 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__SCALE3))
    # --- MTRIX1 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__MTRIX1))
    # --- MTRIX2 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__MTRIX2))
    # --- MTRIX3 ---
    write_records(io, pdb_info, Val(RECORD_TYPE__MTRIX3))
    # --- TVECT ---
    write_records(io, pdb_info, Val(RECORD_TYPE__TVECT))
end

function write_coordinate_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    # first, figure out if we have to write multiple models...
    models = frame_ids(ac)
    
    if pdb_info.selected_model != -1
        models = collect(Iterators.filter(m -> m == pdb_info.selected_model, models))
    end

    multiple_models = length(models) > 1

    for model in models
        if multiple_models
            # write a models record
            write_record(io, pdb_info, RECORD_TAG_MODEL, pdb_info.selected_model)
        end

        # we iterate over all atoms and decide if we have to write an ATOM, HETATM, or TER record

        last_atom = nothing
        atom_number = 1

        for atom in atoms(ac)
            # if the last atom was in a different chain, we have to write a TER record
            if !isnothing(last_atom) && (parent_chain(last_atom) !== parent_chain(atom))
                write_record(io, pdb_info, RECORD_TAG_TER, atom_number, atom_details(last_atom)...)        
                atom_number += 1
            end 

            tag = (isnothing(parent_fragment(atom)) || is_hetero_atom(atom)) ? RECORD_TAG_HETATM : RECORD_TAG_ATOM

            # normalize the name
		    #  if the atom name starts with the element name and the element
		    #  name is a single character (e.g. C, N, O, H and the name is not
		    #  prefixed by a number) then the name should pe prefixed by
		    #  a blank to distinguish CA (carbon alpha) from CA (calcium)
            atom_name = strip(atom.name)
            element_symbol = string(atom.element)
            alt_loc_id = get_property(atom, :alternate_location_id, "")
        
            if (length(atom_name) < 4) && startswith(atom_name, element_symbol) && (length(element_symbol) == 1)
                atom_name = " " * atom_name
            end
            
            write_record(io, pdb_info, tag, 
                  atom_number, atom_name, alt_loc_id, 
                  atom_details(atom)..., atom.r...,
                  get_property(atom, :occupancy, 1.0),
                  get_property(atom, :tempfactor, 0.0),
                  "", element_symbol, atom.formal_charge > 0 ? string(atom.formal_charge) : "")

            atom_number += 1
            last_atom = atom
        end

        # if there were any atoms, write the final TER record
        if natoms(ac) > 0
            write_record(io, pdb_info, RECORD_TAG_TER, atom_number, atom_details(atoms(ac)[end])...)
        end

        if multiple_models
            # write a model record
            write_record(io, pdb_info, RECORD_TAG_ENDMDL)
        end
    end
end

function write_connectivity_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <:Real}
    function needs_connect(b::Bond)
        if ( has_flag(b, :TYPE__COVALENT) ||
             has_flag(b, :TYPE__UNKNOWN)  ||
             has_flag(b, :TYPE__PEPTIDE)  ||
             has_flag(b, :TYPE__DISULPHIDE_BOND))
        
            a, b = get_partners(b)

            if is_hetero_atom(a) || is_hetero_atom(b)
                return true
            end

            ar = parent_fragment(a)
            br = parent_fragment(b)

            if ( a.element === Elements.S && 
                 b.element === Elements.S &&
                 !isnothing(ar) && 
                 !isnothing(br) &&
                 ar !== br)
               
                # we found a disulphide bond
                # TODO: we should check that we build an SSBOND record, too
                return true
            end

            return false
        elseif has_flag(b, :TYPE__HYDROGEN) || has_flag(b, :TYPE__SALT_BRIDGE)
            return true
        end

        return false
    end

    to_connect = Iterators.filter(needs_connect, bonds(ac))

    # we need to divide the bonds up into groups of 3, all sharing the first atom
    bond_groups = groupby(sort(DataFrame(map(b -> (b.a1, b.a2), to_connect)), :1), :1)

    for bg in bond_groups
        for current_partners in Iterators.partition(bg[!, :2], 3)
            partner_strings = vcat(map(string, current_partners), repeat([""], 3 - length(current_partners)))

            # write a CONECT record
            write_record(io, pdb_info, Val(RECORD_TYPE__CONECT), string(bg[1, :1]), partner_strings...)
        end
    end
end

function write_bookkeeping_section(io::IO, pdb_info::PDBInfo, ac::AbstractAtomContainer{T}) where {T <: Real}
    writer_stats = pdb_info.writer_stats

    write_record(io, pdb_info, RECORD_TAG_MASTER, 
        writer_stats.remark_records, 0,
        writer_stats.het_records, writer_stats.helix_records,
        writer_stats.sheet_records, writer_stats.turn_records,
        writer_stats.site_records, writer_stats.coordinate_transformation_records,
        writer_stats.atomic_coordinate_records, writer_stats.ter_records,
        writer_stats.conect_records, writer_stats.seqres_records
    )
end

function write_pdb(io::IO, ac::AbstractAtomContainer{T}) where {T <:Real}
    # get the PDBInfo object
    pdb_info = get_property(ac, :PDBInfo, PDBInfo{T}())

    pdb_info.writer_stats = PDBWriterStats()

    write_title_section(io, pdb_info)
    write_primary_structure_section(io, pdb_info, ac)
    write_heterogen_section(io, pdb_info)
    write_secondary_structure_section(io, pdb_info, ac)
    write_connectivity_annotation_section(io, pdb_info, ac)
    write_miscellaneous_features_section(io, pdb_info)
    write_crystallographic_section(io, pdb_info)
    write_coordinate_section(io, pdb_info, ac)
    write_connectivity_section(io, pdb_info, ac)
    write_bookkeeping_section(io, pdb_info, ac)
    write_record(io, pdb_info, RECORD_TAG_END)
end