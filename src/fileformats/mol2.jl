using BiochemicalAlgorithms: Molecule, Atom, BondOrder, BondOrderType, Bond, Elements, Element, Vector3

 export load_mol2, export_mol2
 
 # TODO: functions for import and export of mol2 files 

 function export_mol2()
    
 end


 function load_mol2(fname::AbstractString, T = Float32)
    mol2_file = open(fname)
    new_mol = Molecule(fname)
    section = ""
    section_line = 0

    for (i, line) in enumerate(readlines(mol2_file))
        if !isempty(line)
            if line[1] == '@'
                println(line)
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
            number = parse(Int64, line[1:7])
            name = strip(line[9:15])
            element = parse(Elements, mol2_get_element(name))
            atomtype = strip(line[51:54])
            r = Vector3{T}(parse(T, line[17:27]), parse(T, line[28:38]), parse(T, line[39:49]))
            v = Vector3{T}(0.0, 0.0, 0.0)
            F = Vector3{T}(0.0, 0.0, 0.0)
            has_velocity = false
            has_force = false
            frame_id = 1
            new_atom = (number = number, name = name, element = element, atomtype = atomtype, r = r, 
                        v = v, F = F, has_velocity = has_velocity, has_force = has_force, frame_id = frame_id)
            push!(new_mol.atoms, new_atom)
        elseif section == "@<TRIPOS>BOND"
            section_line += 1
            new_bond = (a1 = parse(Int64, line[8:13]),
                        a2 = parse(Int64, line[14:19]),
                        order = BondOrderType(mol2_get_BondOrder(line)))
            push!(new_mol.bonds, new_bond)
        end
    end
    return new_mol
 end


 function mol2_get_BondOrder(line::AbstractString)
    bondorder_string = strip(line[20:lastindex(line)])
    if isnumeric(bondorder_string[1])
        return Int(DefBond(parse(Int64, bondorder_string)))
    else
        return Int(parse(DefBonds, bondorder_string))
    end
 end


 function mol2_get_element(name::AbstractString)
    common_elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    first_letter = string(name[1])
    two_letters = ""
    if lastindex(name) >= 2 && !isnumeric(name[2]) 
        two_letters = string(name[1], lowercase(name[2])) 
    end
    if in(common_elements).(two_letters)
        element = two_letters
    elseif in(common_elements).(first_letter)
        element = first_letter
    end
    return element
 end