using Scanf

using BiochemicalAlgorithms

function parse_element_string(es, atom_name)
    result = Elements.Unknown

    # handle special cases
    if es == "D"
        result = Elements.H
    elseif es == "X"
        result = Elements.Unknown
    else
        result = tryparse(ElementType, es)
        
        # if the element symbol is present and valid, it has precedence,
		# otherwise, we try to extract it from the atom name
        # NOTE: this is dangerous if non-PDB names are present; names such
        #       as HE12 are interpreted as Elements.He!
        if isnothing(result)
            symbol =
                if atom_name[1] == ' ' || isdigit(atom_name[1])
                    if atom_name[2] == ' '
                        ""
                    else
                        atom_name[2:2]
                    end 
                else
                    atom_name[1:2]
                end

            result = tryparse(ElementType, symbol)

            if isnothing(result)
                result = Elements.Unknown
                @warn "BiochemicalAlgorithms::PDB::parse_element_string: could not parse element from $(es), $(atom_name); returning Unknown"
            end
        end
    end

    return result
end

function handle_record(line, sys, pdb_info; strict_line_checking, kwargs...)
    # The PDB format description says: "Each line in the PDB entry file
    # consists of 80 columns."
    if length(line) ≠ 80
        if strict_line_checking
            @error "Invalid record: line has $(length(line)) characters instead of 80! Try reading with strict_line_checking=false"
        elseif length(line) < 80
            line *= repeat(" ", 80-length(line))
        end
    end

    # find the record type corresponding to this line
    tag = line[1:6]

    if tag ∉ keys(RECORD_MAP)
        # handle unknown record type
        @warn "Record of PDB not identified: $tag. Record will be skipped." 
    else
        parse_record(
            line, 
            RECORD_MAP[tag]; 
            sys=sys, 
            pdb_info=pdb_info, 
            kwargs...)
    end
end

function ftype(s, T=Float32)
    result =
        if s == "s"
            String
        elseif s == "ld" || s == "d"
            Int
        elseif s == "f"
            T
        elseif s == "c"
            String
        else
            Nothing 
        end
    return result
end

function parse_record(line, record_type; sys::System{T}, pdb_info::PDBInfo{T}, kwargs...) where {T}
    fstring = record_type.format_string
    results = []

    # we start after the tag, which is 6 letters long
    group_start = 7
    for m in eachmatch(r"(\s*)%(-?)(\d+)?(\.\d+)?([[:alpha:]]+)", fstring)
        t = ftype(m.captures[5], T)

        # skip preceeding whitespace
        group_start += length(m.captures[1])

        group_end = group_start +
                    (isnothing(m.captures[3]) ? 1 : parse(Int, m.captures[3]))

        s = line[group_start:group_end-1]

        group_start = group_end

        if t == String
            append!(results, [s])
        else
            _, v = Scanf.scanf(s, Scanf.Format("%" * m.captures[5]), t)
            append!(results, [v])
        end
    end

    interpret_record(Val(record_type.record_type), record_type.record_type, results...; sys=sys, pdb_info, kwargs...)
end

function interpret_record(::Val{RECORD_TYPE__HEADER}, tag, classification, deposition_date, id; sys, pdb_info, kwargs...)
    pdb_info.name = classification
    pdb_info.deposition_date = deposition_date
    pdb_info.id = id
end

function interpret_record(::Val{RECORD_TYPE__TITLE}, tag, continuation, title; sys, pdb_info, kwargs...)
    pdb_info.title *= title
end

function interpret_record(::Val{RECORD_TYPE__TER}, tag, serial_number, residue_name, chain_id,
    residue_number, residue_insertion_code; sys, pdb_info, kwargs...)

    # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # close the current chain
    pdb_info.current_chain = nothing
end

function interpret_record(
    ::Val{RECORD_TYPE__ATOM},
    tag,
    serial_number,
    atom_name,
    alternate_location_identifier,
    residue_name,
    chain_id,
    residue_sequence_number,
    residue_insertion_code,
    x,
    y,
    z,
    occupancy,
    temperature_factor,
    segment_id,
    element_symbol,
    charge;
    sys::System{T},
    pdb_info::PDBInfo{T},
    ignore_xplor_pseudo_atoms,
    is_hetero_atom=false,
    kwargs...) where {T}

    # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # is this an XPLOR pseudo atom that we should skip?
    if ignore_xplor_pseudo_atoms && x >= 9998.0 && y >= 9998.0 && z >= 9998.0
        return
    end

    # is this a new chain?
    if (isnothing(pdb_info.current_chain)
        ||
        chain_id != pdb_info.current_chain.name)
        pdb_info.current_chain = Chain(molecules(sys)[1]; name=chain_id)
        pdb_info.current_residue = nothing
    end

    # right now, we only read the first alternate location indicator
    if (strip(alternate_location_identifier) != ""
        &&
        (alternate_location_identifier != pdb_info.alternate_location_identifier))
        return
    end

    # is this a new residue?
    if (isnothing(pdb_info.current_residue)
        || residue_sequence_number != pdb_info.current_residue.number
        || residue_name != pdb_info.current_residue.name
        || residue_insertion_code != get(pdb_info.current_residue.properties, :insertion_code, ""))

        pdb_info.current_residue = Fragment(
            pdb_info.current_chain,
            residue_sequence_number;
            name=residue_name,
            variant=is_hetero_atom ? FragmentVariant.None : FragmentVariant.Residue,
            properties=Properties([
                :is_hetero_fragment => is_hetero_atom,
                :insertion_code => residue_insertion_code,
                :alternate_location_id => alternate_location_identifier
            ])
        )
    end

    # finally, create the atom
    formal_charge = tryparse(Int, charge)
    formal_charge = isnothing(formal_charge) ? 0 : formal_charge

    Atom(
        pdb_info.current_residue,
        serial_number,
        parse_element_string(strip(element_symbol), atom_name);
        name=String(strip(atom_name)),
        r=Vector3{T}(x, y, z),
        formal_charge=formal_charge,
        properties=Properties([
            :tempfactor => temperature_factor,
            :occupancy => occupancy,
            :is_hetero_atom => is_hetero_atom,
            :insertion_code => residue_insertion_code,
            :alternate_location_id => alternate_location_identifier
        ]),
        flags=Flags(),
        frame_id=pdb_info.current_model
    )
end


function interpret_record(
    ::Val{RECORD_TYPE__HETATM},
    tag,
    serial_number,
    atom_name,
    alternate_location_identifier,
    residue_name,
    chain_id,
    residue_sequence_number,
    residue_insertion_code,
    x,
    y,
    z,
    occupancy,
    temperature_factor,
    segment_id,
    element_symbol,
    charge;
    sys::System{T},
    pdb_info::PDBInfo{T},
    kwargs...) where {T}

    # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # hand over parsing to ATOM
    interpret_record(Val(RECORD_TYPE__ATOM), tag, serial_number,
        atom_name, alternate_location_identifier, residue_name,
        chain_id, residue_sequence_number, residue_insertion_code,
        x, y, z, occupancy, temperature_factor, segment_id,
        element_symbol, charge;
        sys=sys, pdb_info=pdb_info, is_hetero_atom=true, kwargs...)

    # ensure the residue is marked as a hetero fragement, even if the
    # first atom was a regular one
    set_property!(pdb_info.current_residue, :is_hetero_fragment, true)
    set_property!(pdb_info.current_residue, :is_non_standard, true)
    pdb_info.current_residue.variant = FragmentVariant.None

    # recognize water
    if !isnothing(match(r"^OHH|HOH|HHO|H2O|2HO|OH2|SOL|TIP|TIP2|TIP3|TIP4|WAT|D2O$", pdb_info.current_residue.name))
        set_property!(pdb_info.current_residue, :is_water, true)
    end

end

function interpret_record(
    ::Val{RECORD_TYPE__SSBOND},
    tag,
    serial_number,
    first_partner_name,
    first_partner_chain_id,
    first_partner_number,
    first_partner_insertion_code,
    second_partner_name,
    second_partner_chain_id,
    second_partner_number,
    second_partner_insertion_code,
    symmetry_operator_0,
    symmetry_operator_1;
    sys,
    pdb_info,
    kwargs...)

    push!(pdb_info.ssbonds, SSBondRecord(serial_number,
        UniqueResidueID(
            first_partner_name,
            first_partner_chain_id,
            first_partner_number,
            first_partner_insertion_code
        ),
        UniqueResidueID(
            second_partner_name,
            second_partner_chain_id,
            second_partner_number,
            second_partner_insertion_code
        )))
end

function interpret_record(
    ::Val{RECORD_TYPE__HELIX},
    tag,
    serial_number,
    helix_id, 
    initial_residue_name,
    initial_residue_chain_id,
    initial_residue_number,
    initial_residue_insertion_code, 
    terminal_residue_name,
    terminal_residue_chain_id,
    terminal_residue_number,
    terminal_residue_insertion_code,
    helix_class,
    comment,
    length;
    sys,
    pdb_info,
    kwargs...)

    push!(pdb_info.secondary_structures, HelixRecord(
        serial_number,
        helix_id,
        UniqueResidueID(
            initial_residue_name,
            initial_residue_chain_id,
            initial_residue_number,
            initial_residue_insertion_code
        ),
        UniqueResidueID(
            terminal_residue_name,
            terminal_residue_chain_id,
            terminal_residue_number,
            terminal_residue_insertion_code
        ),
        helix_class,
        comment)
    )
end

function interpret_record(
    ::Val{RECORD_TYPE__SHEET},
    tag,
    serial_number,
    sheet_id,
    number_of_strands,
    initial_residue_name,
    initial_residue_chain_id,
    initial_residue_number,
    initial_residue_insertion_code,
    terminal_residue_name,
    terminal_residue_chain_id,
    terminal_residue_number,
    terminal_residue_insertion_code,
    sense_of_strand,
    atom_name_in_current_strand,
    residue_in_current_strand_name,
    residue_in_current_strand_chain_id,
    residue_in_current_strand_number,
    residue_in_current_strand_insertion_code,
    atom_name_in_previous_strand,
    residue_in_previous_strand_name,
    residue_in_previous_strand_chain_id,
    residue_in_previous_strand_number,
    residue_in_previous_strand_insertion_code;
    sys,
    pdb_info,
    kwargs...)

    push!(pdb_info.secondary_structures, SheetRecord(
        serial_number,
        sheet_id,
        UniqueResidueID(
            initial_residue_name,
            initial_residue_chain_id,
            initial_residue_number,
            initial_residue_insertion_code
        ),
        UniqueResidueID(
            terminal_residue_name,
            terminal_residue_chain_id,
            terminal_residue_number,
            terminal_residue_insertion_code
        ),
        sense_of_strand != 0)
    )
end

function interpret_record(
    ::Val{RECORD_TYPE__TURN},
    tag,
    serial_number,
    turn_id,
    initial_residue_name,
    initial_residue_chain_id,
    initial_residue_number,
    initial_residue_insertion_code,
    terminal_residue_name,
    terminal_residue_chain_id,
    terminal_residue_number,
    terminal_residue_insertion_code,
    comment;
    sys,
    pdb_info,
    kwargs...)

    push!(pdb_info.secondary_structures, TurnRecord(
        serial_number,
        turn_id,
        UniqueResidueID(
            initial_residue_name,
            initial_residue_chain_id,
            initial_residue_number,
            initial_residue_insertion_code
        ),
        UniqueResidueID(
            terminal_residue_name,
            terminal_residue_chain_id,
            terminal_residue_number,
            terminal_residue_insertion_code
        ),
        comment)
    )
end


"""
    _atom_by_number(
        ac::AbstractAtomContainer{T} = default_system(),
        idx::Int
    ) -> Atom{T}

Returns the first `Atom{T}` associated with the given `number` in `sys`. Throws a `KeyError` if no such
atom exists.
"""
@inline function _atom_by_number(
    ac::AbstractAtomContainer{T},
    serial_number::Int;
) where T
    idx = filter(atom -> atom.number == serial_number, atoms(ac)).idx
    isempty(idx) ? nothing : atom_by_idx(parent(ac), first(idx))
end

function interpret_record(
    ::Val{RECORD_TYPE__CONECT},
    tag,
    serial_number,
    bond_atom1,
    bond_atom2,
    bond_atom3,
    bond_atom4,
    hbond_atom1,
    hbond_atom2, 
    salt_bridge_atom1, 
    hbond_atom3, 
    hbond_atom4,
    salt_bridge_atom2;
    sys,
    pdb_info,
    kwargs...)

    # find the corresponding atom which bonds we want to recreate
    nonmetals = [
        Elements.Sb, Elements.As, Elements.At, Elements.Fm, Elements.Ge,
        Elements.H, Elements.Ne, Elements.O, Elements.P, Elements.Po,
        Elements.Rn, Elements.Si, Elements.Te, Elements.Ar, Elements.B, 
        Elements.Br, Elements.C, Elements.Cl, Elements.F, Elements.He, 
        Elements.I, Elements.Kr, Elements.N, Elements.S, Elements.Xe]

 
    a = _atom_by_number(sys, serial_number)
    bond_atoms = [bond_atom1, bond_atom2, bond_atom3]
    hbond_atoms = [hbond_atom1, hbond_atom2, hbond_atom3, hbond_atom4]
    salt_bridge_atoms = [salt_bridge_atom1, salt_bridge_atom2]

    if a.element in nonmetals
        for b_idx in bond_atoms
            b = _atom_by_number(sys, b_idx)
            if !isnothing(b) && !is_bound_to(b, a)
                if b.element in nonmetals
                    flags = Flags()
                    push!(flags, :TYPE__COVALENT)
                    Bond(sys, a.idx, b.idx, BondOrder.Single; flags)
                end
            end
        end
        for h_idx in hbond_atoms
            h = _atom_by_number(sys, h_idx)
            if !isnothing(h) && h.element in nonmetals
                if !is_bound_to(h, a)
                    flags = Flags()
                    push!(flags, :TYPE__HYDROGEN)
                    Bond(sys, a.idx, h.idx, BondOrder.Single; flags)
                end
            end
        end
    end

    # create salt bridges
    for b_idx in salt_bridge_atoms
        b = _atom_by_number(sys, b_idx)
        if !isnothing(b) && !is_bound_to(b, a)
            flags = Flags()
            push!(flags, :TYPE__SALT_BRIDGE)
            Bond(sys, a.idx, b.idx, BondOrder.Single; flags)
        end
    end
end


function interpret_record(record_type, tag, record_data...; sys, pdb_info, kwargs...)
    push!(pdb_info.records, PDBRecord(tag, record_data))
end

function postprocess_ssbonds_!(sys, pdb_info, fragment_cache)
    for ssbond in pdb_info.ssbonds
        f1 = fragment_cache[ssbond.first]
        f2 = fragment_cache[ssbond.second]

        set_flag!(f1, :HAS_SSBOND)
        set_flag!(f2, :HAS_SSBOND)

        # build the bond
        a1 = findfirst(a -> a.element == Elements.S, atoms(f1))
        a2 = findfirst(a -> a.element == Elements.S, atoms(f2))
     
        if !isnothing(a1) && !isnothing(a2)
            # if the information about disulfid bonds were in connections, the bond already exists
            # we will only add the corresponding flag 
            bond_idx = findfirst(
                b -> 
                ((b.a1 == atoms(f1)[a1].idx) && (b.a2 == atoms(f2)[a2].idx)) ||
                ((b.a1 == atoms(f2)[a2].idx) && (b.a2 == atoms(f1)[a1].idx)), bonds(sys))
           
            if !isnothing(bond_idx)
                set_flag!(bonds(sys)[bond_idx],(:TYPE__DISULPHIDE_BOND)) 
            else
                Bond(sys, atoms(f1)[a1].idx, atoms(f2)[a2].idx, BondOrder.Single; flags=Flags((:TYPE__DISULPHIDE_BOND,)))
            end
        end
    end
end

function postprocess_secondary_structures_!(sys, pdb_info, fragment_cache, create_coils)
    for ss in pdb_info.secondary_structures
        initial_res  = fragment_cache[ss.initial_residue]
        terminal_res = fragment_cache[ss.terminal_residue]

        if parent_chain(initial_res) != parent_chain(terminal_res)
            @error "Invalid Secondary Structure record: residues in different chains!"
        end

        new_ss = if typeof(ss) == HelixRecord
            new_ss = SecondaryStructure(
                parent_chain(initial_res), 
                ss.number, 
                SecondaryStructureElement.Helix;
                name=ss.name)
            set_property!(new_ss, :HELIX_CLASS, ss.helix_class)
            set_property!(new_ss, :COMMENT, ss.comment)

            new_ss
        elseif typeof(ss) == SheetRecord
            new_ss = SecondaryStructure(
                parent_chain(initial_res), 
                ss.number, 
                SecondaryStructureElement.Strand;
                name="$(ss.name):$(ss.number)")
            set_property!(new_ss, :STRAND_SENSE, ss.sense_of_strand)

            new_ss
        else
            new_ss = SecondaryStructure(
                parent_chain(initial_res), 
                ss.number, 
                SecondaryStructureElement.Turn;
                name=ss.name)
            set_property!(new_ss, :COMMENT, ss.comment)

            new_ss
        end

        for f in fragments(parent_chain(initial_res))
            if f.number >= initial_res.number && f.number <= terminal_res.number
                f.secondary_structure_idx = new_ss.idx
            end
        end
    end

    # do we need to put every amino acid residue into a secondary structure?
    if create_coils
        for c in chains(sys)
            fs = fragments(c)

            if isempty(fs[fs.secondary_structure_idx .== nothing])
                continue
            end

            last_fragment = nothing

            current_ss = nothing
            for f in fs
                if isnothing(f.secondary_structure_idx)
                    if isnothing(last_fragment) || parent_secondary_structure(last_fragment) != current_ss
                        current_ss = SecondaryStructure(c, 
                        maximum(secondary_structures(c).number, init=0) + 1, 
                        SecondaryStructureElement.Coil)
                    end
                    f.secondary_structure_idx = current_ss.idx
                end
                last_fragment = f
            end
        end
    end

    # finally, renumber the elements to be consecutive    
    sort_secondary_structures!(sys, by=se -> (se.chain_idx, minimum(f.number for f in fragments(se))))
    
    for c in chains(sys)
        secondary_structures(c).number .= collect(1:nsecondary_structures(c))
    end
end

function load_pdb(
        filename::String, 
        T;
        keep_metadata=true,
        strict_line_checking=true,
        selected_model=-1, 
        ignore_xplor_pseudo_atoms=true,
        create_coils=true)
    pdblines = readlines(filename)

    sys = System{T}("")
    mol = Molecule(sys; name="")

    pdb_info = PDBInfo{T}(selected_model)

    for pl in pdblines
        handle_record(pl, sys, pdb_info; 
            strict_line_checking=strict_line_checking,
            ignore_xplor_pseudo_atoms=ignore_xplor_pseudo_atoms)
    end

    sys.name = strip(pdb_info.name)
    mol.name = sys.name

    fragment_cache = Dict{UniqueResidueID, Fragment{T}}()

    for f in fragments(sys)
        fragment_cache[UniqueResidueID(f.name, parent_chain(f).name, f.number, get_property(f, :insertion_code, " "))] = f
    end

    postprocess_ssbonds_!(sys, pdb_info, fragment_cache)
    postprocess_secondary_structures_!(sys, pdb_info, fragment_cache, create_coils)

    if keep_metadata
        set_property!(sys, :PDBInfo, pdb_info)
    end

    sys
end
