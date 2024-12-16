export
    load_hinfile,
    write_hinfile

mutable struct HINParserState{T}
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

function _handle_hin_atom!(parser_state::HINParserState{T}) where T
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

function _handle_hin_velocity!(parser_state::HINParserState{T}) where T
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

function _handle_hin_residue!(parser_state::HINParserState)
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

function _handle_hin_endresidue!(parser_state::HINParserState)
    if parser_state.state ≠ :IN_RESIDUE
        @error "illegal HIN file in line $(parser_state.line_number): <endres> tag may appear only inside a <res> tag!"
    end

    parser_state.state = :IN_MOLECULE

    parser_state.current_fragment = nothing
end

function _handle_hin_molecule!(parser_state::HINParserState)
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

function _handle_hin_endmolecule!(parser_state::HINParserState)
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

function _handle_hin_system!(parser_state::HINParserState{T}) where T
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

function _handle_hin_box!(parser_state::HINParserState{T}) where T
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

function _handle_hin_line!(parser_state::HINParserState)
    # the HyperChem tag is always the first word in the line
    tag = first(split(parser_state.line))

    if tag == "atom"
        _handle_hin_atom!(parser_state)
    elseif tag == "vel"
        _handle_hin_velocity!(parser_state)
    elseif tag == "res"
        _handle_hin_residue!(parser_state)
    elseif tag == "endres"
        _handle_hin_endresidue!(parser_state)
    elseif tag == "mol"
        _handle_hin_molecule!(parser_state)
    elseif tag == "endmol"
        _handle_hin_endmolecule!(parser_state)
    elseif tag == "sys"
        _handle_hin_system!(parser_state)
    elseif tag == "box"
        _handle_hin_box!(parser_state)
    elseif tag  ∉ ("forcefield", "user1color", "user2color", "user3color", "user4color", "view", "seed",
                   "mass", "basisset", "selection", "endselection", "selectrestraint", "selectatom", "formalcharge") # we know these tags, but ignore them
        @warn "unexpected HIN file format in line $(parser_state.line_number): unkown tag $(tag)!"
    end
end

"""
    load_hinfile(
        fname::String,
        ::Type{T} = Float32
    ) -> System{T}

Read a HyperChem HIN file.
"""
function load_hinfile(fname::String, ::Type{T} = Float32) where {T <: Real}
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

        _handle_hin_line!(parser_state)
    end

    parser_state.s
end

"""
    write_hinfile(
        fname::String,
        ac::AtomContainer
    )

Save an atom container as HyperChem HIN file.

!!! note
    HIN files define molecules as connected components in the molecular graph. If the AtomContainer is
    missing bonds, e.g., after reading a PDB file and not postprocessing it correctly, the HIN file may
    contain a surprisingly large number of molecules.
"""
function write_hinfile(fname::String, ac::AbstractAtomContainer)
    mg = convert(MolecularGraph.SDFMolGraph, ac)
    hin_molecules = connected_components(mg)

    open(fname, "w") do io
        write(io, "; HyperChem file created by BiochemicalAlgorithms\n")
        write(io, ";\n")
        write(io, "forcefield AMBER\n")

        write(io, "sys $(get_property(ac, :temperature, 0.0))\n")

        if has_property(ac, :periodic_box_width) && has_property(ac, :periodic_box_height) && has_property(ac, :periodic_box_depth)
            write(io, "box $(get_property(ac, :periodic_box_width)) $(get_property(ac, :periodic_box_height)) $(get_property(ac, :periodic_box_depth))\n")
        end

        # now, for each connected component, create a molecule
        for (mi, m) in enumerate(hin_molecules)
            # getting the name of the molecule is not entirely simple...
            # we use the first atom of the connected component to map back to the parent molecule, as many components can arise from the same parent molecule
            mol_name = parent_molecule(atom_by_idx(ac, mg.gprops[:atom_idx][first(m)])).name

            # names cannot contain double quotes
            mol_name = replace(strip(mol_name), "\"" => "")

            # newer HyperChem releases require molecule names to be enclosed in double quotes
            write(io, mol_name != "" ? "mol $(mi) \"$(mol_name)\"\n" : "mol $(mi)\n")

            # we need to re-number the atoms; the input are global numbers, unique over all molecules, while HIN wants atom numbers to be relative to the
            # current molecule
            index_map = Dict(a => ai for (ai, a) in enumerate(m))

            # now, handle each atom
            last_fragment = nothing
            num_fragments = 1
            for a in m
                orig_atom = atom_by_idx(ac, mg.gprops[:atom_idx][a])

                current_fragment = parent_fragment(orig_atom)

                if current_fragment != last_fragment
                    # do we need to close a fragment?
                    if !isnothing(last_fragment)
                        write(io, "endres $(num_fragments)\n")
                        num_fragments += 1
                    end

                    # do we need to open a new fragment?
                    if !isnothing(current_fragment)
                        frag_name = strip(current_fragment.name) != "" ? strip(current_fragment.name) : "-"
                        chain_name = strip(parent_chain(current_fragment).name) != "" ? strip(parent_chain(current_fragment).name) : "-"

                        het_id = (is_amino_acid(current_fragment) || current_fragment.variant == FragmentVariant.Residue) ? "-" : "h"
                        write(io, "res $(num_fragments) $(frag_name) $(current_fragment.number) $(het_id) $(chain_name)\n")
                    end

                    last_fragment = current_fragment
                end

                # now, write the atom itself
                atom_name = strip(orig_atom.name) == "" ? "-" : strip(orig_atom.name)
                if occursin(" ", atom_name)
                    @warn "Atom names in HIN files cannot contain spaces!"

                    atom_name = first(split(atom_name))
                end

                atom_type = strip(orig_atom.atom_type) == "" ? "-" : strip(orig_atom.atom_type)

                atom_string  = "atom $(index_map[a]) $(atom_name) $(orig_atom.element) $(atom_type) - $(orig_atom.charge) $(orig_atom.r[1]) $(orig_atom.r[2]) $(orig_atom.r[3])"
                atom_string *= " $(nbonds(orig_atom))"

                for b in neighbors(mg, a)
                    # find the corresponding edge
                    e = mg.eprops[MolecularGraph.u_edge(mg, a, b)]

                    order =
                        if e.order == 1
                            's'
                        elseif e.order == 2
                            'd'
                        elseif e.order == 3
                            't'
                        else
                            '-'
                        end

                    atom_string *= " $(index_map[b]) $(order)"
                end

                atom_string *= "\n"

                write(io, atom_string)

                # write the velocities
                write(io, "vel $(index_map[a]) $(orig_atom.v[1]) $(orig_atom.v[2]) $(orig_atom.v[3])\n")
            end

            write(io, "endmol $(mi)\n")
        end
    end
end
