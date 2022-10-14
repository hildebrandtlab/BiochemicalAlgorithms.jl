
 export load_mol2, export_mol2
 
 # TODO: functions for import and export of mol2 files 

 function export_mol2()
    
 end


 function load_mol2(fname::AbstractString)
    mol2_file = open(fname)
    new_mol = Molecule(fname)
    section = ""
    section_line = 0

    for (i, line) in enumerate(readlines(mol2_file))
        if line[1] == "@"
            section = string(line)
            section_line = 0
            continue
        end
        if section == "@<TRIPOS>MOLECULE"
            section_line += 1
            if section_line == 1
                new_mol.name == string(line)
            end
        elseif section == "@<TRIPOS>ATOM"
            section_line += 1
            new_atom = Atom()
            new_atom.number = parse(Int64, line[1:7])
            new_atom.name = strip(line[9:15])
            new_atom.element = mol2_get_element(new_atom.name)
            new_atom.r = mol2_get_r_vector(line)
            new_atom.v = Ref(Vector3{T}(0.0, 0.0, 0.0))
            new_atom.F = Ref(Vector3{T}(0.0, 0.0, 0.0))
            new_atom.has_velocity = false
            new_atom.has_force = false
            new_atom.frame_id = 1
        elseif section == "@<TRIPOS>BOND"
            section
        end
    end
 end

 function mol2_get_element(name::AbstractString)
    common_elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
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

 function mol2_get_r_vector(line)
    r_vector = Vector3{Float32}
    r_vector = [parse(Float32, line[17:27]), parse(Float32, line[28:38]), parse(Float32, line[39:49])]
    return r_vector
 end