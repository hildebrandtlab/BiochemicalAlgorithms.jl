
using BiochemicalAlgorithms: Molecule, Bond, Atom
using DataFrames, DelimitedFiles, CSV

export Atomtype, make_dict_from_dat

function atomtype_gaff(num::UInt128, mol::Molecule)
    atom_gaff = mol.atoms[num]
    bonds_gaff = subset(mol.bonds, aid1 == num)
    
    for i = 1:length(max(Molecule.bonds.a1))
    end
end

function remove_space_add_comma(line::AbstractString)
    out_str = ""
    i = 1
	while i <= lastindex(line)
        if i > 1 && line[i-1] == ' ' && line[i] != ' ' && line[i] != '&'
                out_str = string(out_str,",")
        end
        if line[i] != ' ' && line[i] != '&' && line[i] != " "
            out_str = string(out_str,line[i])
        end
        if i == lastindex(line)
            out_str = string(out_str,'&')
        end
        i += 1
	end
    return out_str
end

function df_maker(mapfile::AbstractString)
    source_file = open(mapfile)
    temp_str = ""
    temp_arr = []
    for line in readlines(source_file)
        if lastindex(line)>=3 && line[1:3] == "ATD"
            temp_var = remove_space_add_comma(line)
            temp_str = string(temp_str, temp_var)
            append!(temp_arr,temp_str)            
        end
    end
    # temp_str = readdlm(temp_str, ',', AbstractString, '&')
    writedlm("src\\mappings\\antechamber\\atomtype_gaff.csv", temp_arr, ',')
    df_atomtype = CSV.read("src\\mappings\\antechamber\\atomtype_gaff.csv")
    return df_atomtype
end