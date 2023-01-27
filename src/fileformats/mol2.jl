using BiochemicalAlgorithms

export load_mol2, export_mol2


function load_mol2(fname::AbstractString, T = Float32)
    # For validation of small molecules only: no <TRIPOS>SUBSTRUCTURES import implemented
    if fname[lastindex(fname)-4:lastindex(fname)] != ".mol2"
        println("Please make sure, you are loading a mol2 file!\n(file ending should be mol2)")
        return
    end
    new_mol = Molecule(fname)
    section = ""
    section_line = 0

    for (i, line) in enumerate(readlines(fname))
        if !isempty(line)
            if line[1] == '@'
                section = string(line)
                section_line = 0
                continue
            end
        end
        if section == "@<TRIPOS>MOLECULE"
            section_line += 1
            if section_line == 1
                new_mol.name == string(line)
            end
        elseif section == "@<TRIPOS>ATOM"
            section_line += 1
            line_elements = split(line)
            new_atom = (number = parse(Int64, line_elements[1]), 
                        name = line_elements[2], 
                        element = parse(Elements, mol2_get_element(line_elements[2])), 
                        atomtype = line_elements[6], 
                        r = Vector3{T}(parse(T, line_elements[3]), parse(T, line_elements[4]), parse(T, line_elements[5])), 
                        v = Vector3{T}(0.0, 0.0, 0.0), F = Vector3{T}(0.0, 0.0, 0.0), 
                        has_velocity = false, has_force = false, frame_id = parse(Int64, line_elements[7]),
                        properties = Dict{String, Union{T, String}}("Charge" => parse(T, line_elements[9]), "Chain" => string(line_elements[8])))
            push!(new_mol.atoms, new_atom)
        elseif section == "@<TRIPOS>BOND"
            section_line += 1
            line_elements = split(line)
            order_ = BondOrderType(mol2_get_BondOrder(line_elements[4]))
            properties_dict = Dict{String, Union{T, String}}()
            if order_ == BondOrder.Unknown && line_elements[4] == "ar" || 
                order_ == BondOrder.Single && line_elements[4] == "am"
                println(line_elements[4])
                properties_dict["TRIPOS_tag"] = string(line_elements[4])
            end
            new_bond = (a1 = parse(Int64, line_elements[2]),
                        a2 = parse(Int64, line_elements[3]),
                        order = order_, 
                        properties = properties_dict)
            push!(new_mol.bonds, new_bond)
        end
    end
    return new_mol
end


function export_mol2(mol::Molecule, filelocation::String)
    # For validation of small molecules only: no <TRIPOS>SUBSTRUCTURES export implemented
    mol_name = (!isnothing(findlast('.', mol.name)) && in(findlast('.', mol.name), [lastindex(mol.name)-5:lastindex(mol.name);])) ? basename(mol.name[1:findlast('.', mol.name)-1]) : basename(mol.name)
    export_file = open(string(filelocation, mol_name, ".mol2") , "w")
    
    ### Molecule section
    write(export_file, "@<TRIPOS>MOLECULE\n")
    
    molecule_section_line1 = string(mol_name, "\n")

    num_atoms = build_flush_right_string(nrow(mol.atoms), 5)
    num_bonds = build_flush_right_string(nrow(mol.bonds), 6)
    num_subst = build_flush_right_string(0, 6)
    if typeof(mol) == PDBMolecule{Float32} && !isempty(mol.chains[1].fragments)
        num_subst = build_flush_right_string(nrow(mol.chains[1].fragments), 6)
    end
    num_feat = build_flush_right_string(0, 6)
    num_sets = build_flush_right_string(0, 6)
    molecule_section_line2 = string(num_atoms, num_bonds, num_subst, num_feat, num_sets, "\n")

    mol_type = "SMALL" # BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE possible according to Tripos mol2 specification pdf
    if haskey(mol.properties, "Type")
        mol_type = mol.properties["Type"]
    end
    molecule_section_line3 = string(mol_type, "\n")

    charge_type = "NO_CHARGES" # DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN, 
    # GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_ CHARGES, MMFF94_CHARGES, USER_CHARGES
    # possible according to Tripos mol2 specification
    if haskey(mol.properties, "Charge_Type")
        mol_type = mol.properties["Charge_Type"]
    end
    molecule_section_line4 = string(charge_type, "\n")

    status_bits = "\n"
    molecule_section_line5 = status_bits
    mol_comment = "\n"
    molecule_section_line6 = mol_comment

    write(export_file, molecule_section_line1, molecule_section_line2, 
            molecule_section_line3, molecule_section_line4,
            molecule_section_line5, molecule_section_line6) 

    ### Atom section
    write(export_file, "@<TRIPOS>ATOM\n")

    for i = (1:nrow(mol.atoms))
        atom_id = string(build_flush_right_string(mol.atoms.number[i], 7), " ")
        atom_name = build_flush_left_string(mol.atoms.element[i], 6)
        x_coordinate_string = build_Float32_string(mol.atoms.r[i][1], 13, 4)
        y_coordinate_string = build_Float32_string(mol.atoms.r[i][2], 11, 4)
        z_coordinate_string = build_Float32_string(mol.atoms.r[i][3], 11, 4)
        atom_type = string(" ", build_flush_left_string("DU", 6))
        if !isempty(mol.atoms.atomtype[i])
            atom_type = string(" ", build_flush_left_string(mol.atoms.atomtype[i], 6))
        end
        subst_id = string(build_flush_right_string(1, 6), " ")
        subst_name = build_flush_left_string("nan", 6)
        if typeof(mol) == PDBMolecule{Float32}
            if !isempty(mol.atoms.residue_id[i])
                subst_id = string(build_flush_right_string(mol.atoms.residue_id[i], 6), " ")
            end
            if !isempty(mol.atoms.residue_name[i])
                subst_name = build_flush_left_string(mol.atoms.residue_name[i], 6)
            end
        end
        charge = build_Float32_string(0.0, 10, 4)
        if haskey(mol.atoms.properties[i], "Charge")
            charge = build_Float32_string(convert(Float32, mol.atoms.properties[i]["Charge"]), 10, 4)
        end
        # status_bits never set by user, DSPMOD, TYPECOL, CAP, BACKBONE, DICT, ESSENTIAL, 
        # WATER and DIRECT are possible according to Tripos mol2 specification
        atom_section_line = string(atom_id, atom_name, x_coordinate_string, 
                                    y_coordinate_string, z_coordinate_string,
                                    atom_type, subst_id, subst_name, charge, "\n")
        write(export_file, atom_section_line)
    end 

    ### Bond section
    if !isempty(mol.bonds)
        write(export_file, "@<TRIPOS>BOND\n")

        for i = (1:nrow(mol.bonds))
            bond_id = build_flush_right_string(i, 6)
            origin_atom_id = build_flush_right_string(mol.bonds.a1[i], 6)
            target_atom_id = build_flush_right_string(mol.bonds.a2[i], 6)
            bond_type = string(" ", build_flush_left_string(Int(mol.bonds.order[i]), 4))
            # status_bits never set by user, TYPECOL, GROUP, CAP, BACKBONE, DICT and INTERRES 
            # are possible according to Tripos mol2 specification
            bond_section_line = string(bond_id, origin_atom_id, target_atom_id, bond_type, "\n")
            write(export_file, bond_section_line)
        end
    end

    close(export_file)
end


function build_Float32_string(input::AbstractFloat, length::Int64, decimals::Int64)
    col_string = string(input)
    while  (lastindex(col_string)-findfirst('.', col_string)) < decimals
        col_string = string(col_string, "0")
    end
    while lastindex(col_string) < length
        col_string = string(" ", col_string)
    end
    return col_string
end


function build_flush_right_string(input::Any, length::Int64)
   col_string = string(input)
   while lastindex(col_string) < length
       col_string = string(" ", col_string)
   end
   return col_string
end


function build_flush_left_string(input::Any, length::Int64)
    col_string = string(input)
    while lastindex(col_string) < length
        col_string = string(col_string, " ")
    end
    return col_string
end


function mol2_get_BondOrder(substring::AbstractString)
    if isnumeric(substring[1])
        return Int(BondShortOrderType(parse(Int64, substring)))
    elseif isdefined(BondShortOrder, Symbol(substring))
        return Int(getproperty(BondShortOrder, Symbol(substring)))
    elseif substring == "am"
        return Int(BondOrder.Single)
    else
        return Int(BondOrder.Unknown)
    end
end


function mol2_get_element(name::AbstractString)
    common_elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    first_letter = string(name[1])
    two_letters = ""
    if lastindex(name) >= 2 && !isnumeric(name[2]) 
        two_letters = string(name[1], lowercase(name[2])) 
    end
    if in(two_letters, common_elements)
        element = two_letters
    elseif in(first_letter, common_elements)
        element = first_letter
    end
    return element
end