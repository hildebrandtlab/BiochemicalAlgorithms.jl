export load_hinfile

mutable struct HINParserState{T<:Real}
    line::String
    line_number::Int
    s::Union{System{T}, Nothing}

    state::Symbol

    current_molecule::Union{Molecule{T}, Nothing}
    current_fragment::Union{Fragment{T}, Nothing}
    
    current_bonds::Deque{Tuple{Int, Int, BondOrderType}}

    atom_cache::Dict{Int, Atom{T}}

    function HINParserState{T}() where {T<:Real}
        new("", -1, nothing, :START, nothing, nothing, Deque{Tuple{Int, Int, BondOrderType}}(), Dict{Int, Atom{T}}())
    end
end

function handle_hin_atom_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ∉ (:IN_MOLECULE, :IN_RESIDUE)
        @error "illegal HIN file in line $(parser_state.line_number): <atom> tag may appear only inside a <mol> or <res> tag!"
    end

    fields = split(parser_state.line)

    if length(fields) < 11
        @error "illegal HIN file in line $(parser_state.line_number): <atom> line requires at least 10 arguments!"
    end

    atom_number  = try 
        parse(Int, fields[2])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[2]) as atom number!"
    end

    atom_name = String(fields[3])

    # Amber writes "lone pair atoms" with pseudo-element "Lp" -- we treat these differently; right now,
    # we keep them around and remove them later
    atom_element = fields[4] == "Lp" ? Elements.Unknown : parse(ElementType, fields[4])

    # TODO: BALL used to set the radius to the van der Waals - Radius of the Element; should we do that, too?
    
    atom_type = fields[5] == "**" ? "?" : String(strip(fields[5]))

    charge = try
        parse(T, fields[7])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[7]) as a charge!"
    end

    x, y, z = try
        parse(T, fields[8]), parse(T, fields[9]), parse(T, fields[10])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse ($(fields[8]), $(fields[9]), $(fields[10])) as position vector!"
    end

    ac = if parser_state.state == :IN_RESIDUE
        # PDB HETATMS are denoted by putting an 'h' in field 6
        # TODO: we should use these properties somewhere...
        if 'h' ∈ fields[6]
            unset_flag!(parser_state.current_fragment, :IS_AMINO_ACID)
            set_flag!(parser_state.current_fragment, :NON_STANDARD_RESIDUE)
            parser_state.current_fragment.variant = FragmentVariant.None
        else
            set_flag!(parser_state.current_fragment, :IS_AMINO_ACID)
            unset_flag!(parser_state.current_fragment, :NON_STANDARD_RESIDUE)
            parser_state.current_fragment.variant = FragmentVariant.Residue
        end

        parser_state.current_fragment
    else
        parser_state.current_molecule
    end

    at = Atom(ac, atom_number, atom_element; name=atom_name, atom_type=atom_type, r=Vector3{T}(x, y, z), charge=charge)

    parser_state.atom_cache[atom_number] = at

    num_bonds = try
        parse(Int, fields[11])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[11]) as the number of bonds!"
    end

    if num_bonds > 0
        # check whether the number of fields is consistent with the number of bonds
        if length(fields) != (11 + 2 * num_bonds)
            @error "illegal HIN file in line $(parser_state.line_number): inconsistent number of fields ($(length(fields))) for $(num_bonds) bonds!"
        end

        for current_bond in range(1, num_bonds)
            partner_number = try
                parse(Int, fields[11 + 2*current_bond - 1])
            catch _
                @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[11 + 2*current_bond - 1]) as atom index!"
            end

            partner_type_field = strip(String(fields[11 + 2*current_bond])) 
            partner_type = 
                if partner_type_field == 's'
                    BondOrder.Single
                elseif partner_type_field == 'd'
                    BondOrder.Double
                elseif partner_type_field == 't'
                    BondOrder.Triple
                elseif partner_type_field == 'a'
                    BondOrder.Aromatic
                else
                    BondOrder.Unknown
                end

            push!(parser_state.current_bonds, (atom_number, partner_number, partner_type))
        end
    end
end

function handle_hin_velocity_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ∉ (:IN_MOLECULE, :IN_RESIDUE)
        @error "illegal HIN file in line $(parser_state.line_number): <vel> tag may appear only inside a <mol> or <res> tag!"
    end

    fields = split(parser_state.line)

    atom_number  = try 
        parse(Int, fields[2])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[2]) as atom number!"
    end
    
    vel_x, vel_y, vel_z = try
        parse(T, fields[3]), parse(T, fields[4]), parse(T, fields[5])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse ($(fields[3]), $(fields[4]), $(fields[5])) as velocity vector!"
    end

    # the atom must already exist
    if atom_number ∉ keys(parser_state.atom_cache)
        @error "illegal HIN file in line $(parser_state.line_number): cannot store velocity for atom $(atom_number), which has not been seen previously!"
    end

    parser_state.atom_cache[atom_number].v = Vector3{T}(vel_x, vel_y, vel_z)
end

function handle_hin_residue_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :IN_MOLECULE
        @error "illegal HIN file in line $(parser_state.line_number): <res> tag may appear only inside a <mol> tag!"
    end

    fields = split(parser_state.line)

    parser_state.state = :IN_RESIDUE

    # do we have a chain already?
    if nchains(parser_state.current_molecule) == 0
        Chain(parser_state.current_molecule)
    end

    chain = first(chains(parser_state.current_molecule))

    residue_name = strip(fields[3]) == "-" ? "" : String(strip(fields[3]))

    residue_number = try
        parse(Int, fields[4])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[4]) as residue number!"
    end

    chain_name = strip(String(fields[6])) == "-" ? "" : fields[6]

    chain.name = chain_name

    parser_state.current_fragment = Fragment(chain, residue_number; name=residue_name)
end

function handle_hin_endresidue_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :IN_RESIDUE
        @error "illegal HIN file in line $(parser_state.line_number): <endres> tag may appear only inside a <res> tag!"
    end

    parser_state.state = :IN_MOLECULE

    parser_state.current_fragment = nothing
end

function handle_hin_molecule_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :START
        @error "illegal HIN file in line $(parser_state.line_number): <mol> tag may appear only on the top level!"
    end

    fields = split(parser_state.line)

    molecule_name = ""
    if length(fields) > 2
        molecule_name = join(fields[3:end], " ")
    end

    molecule_name = replace(molecule_name, "\"" => "")

    parser_state.state = :IN_MOLECULE

    parser_state.current_molecule = Molecule(parser_state.s; name=molecule_name)
end

function handle_hin_endmolecule_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :IN_MOLECULE
        @error "illegal HIN file in line $(parser_state.line_number): <endmol> tag may appear only inside a <mol> tag!"
    end

    parser_state.state = :START

    # do we have an empty fragment to clean up?
    if !isnothing(parser_state.current_fragment) && natoms(parser_state.current_fragment) == 0
        delete!(parser_state.current_fragment)
    end
    parser_state.current_fragment = nothing

    # or an empty chain?
    if nchains(parser_state.current_molecule) >= 1 && natoms(first(chains(parser_state.current_molecule))) == 0
        delete!(first(chains(parser_state.current_molecule)))
    end

    # now, build all bonds
    for b in parser_state.current_bonds
        a1 = get(parser_state.atom_cache, b[1], nothing)
        a2 = get(parser_state.atom_cache, b[2], nothing)
        
        if nothing ∈ (a1, a2)
            @error "illegal HIN file in line $(parser_state.line_number): cannot finalize bond $(b)!"
        end

        # did we already process this bond?
        if is_bound_to(a1, a2)
            continue
        end

        type = b[3]

        nb = Bond(parser_state.s, a1.idx, a2.idx, type)
    
        f1 = parent_fragment(a1)
        f2 = parent_fragment(a2)

        # fix disulphide bridges
        if    a1.element == Elements.S &&
              a2.element == Elements.S &&
              f1 != f2 &&
              !isnothing(f1) &&
              !isnothing(f2) &&
              has_flag(f1, :IS_AMINO_ACID) &&
              has_flag(f2, :IS_AMINO_ACID)
            
            set_flag!(f1, :HAS_SSBOND)
            set_flag!(f2, :HAS_SSBOND)

            set_flag!(nb, :DISULPHIDE_BOND)
        end
    end

    empty!(parser_state.current_bonds)
    empty!(parser_state.atom_cache)

    parser_state.current_molecule = nothing
end

function handle_hin_system_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :START
        @error "illegal HIN file in line $(parser_state.line_number): <sys> tag may appear only on the top level!"
    end

    fields = split(parser_state.line)

    temperature = try
        parse(T, fields[2])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse $(fields[2]) as temperature!"
    end

    set_property!(parser_state.s, :temperature, temperature)
end

function handle_hin_box_!(parser_state::HINParserState{T}) where {T<:Real}
    if parser_state.state ≠ :START
        @error "illegal HIN file in line $(parser_state.line_number): <box> tag may appear only on the top level!"
    end

    fields = split(parser_state.line)

    box_width, box_height, box_depth = try
        parse(T, fields[2]), parse(T, fields[3]), parse(T, fields[4])
    catch _
        @error "illegal HIN file in line $(parser_state.line_number): could not parse ($(fields[2]), $(fields[3]), $(fields[4])) as velocity vector!"
    end

    set_property!(parser_state.s, :periodic_box_width,  box_width)
    set_property!(parser_state.s, :periodic_box_height, box_height)
    set_property!(parser_state.s, :periodic_box_depth,  box_depth)
end

function handle_hin_line_!(parser_state::HINParserState{T}) where {T<:Real}
    # the HyperChem tag is always the first word in the line
    tag = first(split(parser_state.line))

    if tag == "atom"
        handle_hin_atom_!(parser_state)
    elseif tag == "vel"
        handle_hin_velocity_!(parser_state)
    elseif tag == "res"
        handle_hin_residue_!(parser_state)
    elseif tag == "endres"
        handle_hin_endresidue_!(parser_state)
    elseif tag == "mol"
        handle_hin_molecule_!(parser_state)
    elseif tag == "endmol"
        handle_hin_endmolecule_!(parser_state)
    elseif tag == "sys"
        handle_hin_system_!(parser_state)
    elseif tag == "box"
        handle_hin_box_!(parser_state)
    elseif tag  ∉ ("forcefield", "user1color", "user2color", "user3color", "user4color", "view", "seed",
                   "mass", "basisset", "selection", "endselection", "selectrestraint", "selectatom", "formalcharge") # we know these tags, but ignore them
        @warn "unexpected HIN file format in line $(parser_state.line_number): unkown tag $(tag)!"
    end
end

"""
    load_hinfile(
        fname::String,
        T=Float32
    ) -> System{T}

Read a HyperChem HIN file.
"""
function load_hinfile(fname::String, T=Float32)
    # we use a very simple state machine for parsing HIN files; legal states are
    #   - :START
    #   - :IN_MOLECULE
    #   - :IN_RESIDUE
    # and legal transitions are
    #   - :START       -> :IN_MOLECULE
    #   - :IN_MOLECULE -> :IN_RESIDUE
    #   - :IN_RESIDUE  -> :IN_MOLECULE
    #   - :IN_MOLECULE -> :START

    parser_state = HINParserState{T}()

    parser_state.state = :START
    parser_state.s = System{T}()

    parser_state.current_molecule = nothing
    parser_state.current_fragment = nothing

    for (i, line) in enumerate(readlines(fname))
        # ignore comment lines
        if line == "" || first(line) == ';'
            continue
        end

        parser_state.line = line
        parser_state.line_number = i

        handle_hin_line_!(parser_state)
    end

    parser_state.s
end