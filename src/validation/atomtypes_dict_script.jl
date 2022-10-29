using BiochemicalAlgorithms

export load_all_gaff_paper_files, run_atomtyping_on_gaff_paper_files

function load_all_gaff_paper_files()
    file_location_a = "test/data/gaff_paper_examples/a/"
    file_location_b = "test/data/gaff_paper_examples/b/"
    row_count_a = lastindex(readdir(file_location_a))
    row_count_total = lastindex(readdir(file_location_a)) + lastindex(readdir(file_location_b))
    mol_df = DataFrame([Vector{String}(undef, row_count_total), Vector{AbstractMolecule}(undef, row_count_total)], ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location_a))
        mol_df.molname[num] = string("mol_a_", i[1:2])
        mol_df.abstract_mol[num] = load_pubchem_json(string(file_location_a, i))
    end
    for (num, i) in enumerate(readdir(file_location_b))
        num_b = num + row_count_a
        mol_df.molname[num_b] = string("mol_b_", i[1:2])
        mol_df.abstract_mol[num_b] = load_pubchem_json(string(file_location_b, i))
    end
    return mol_df
end

function run_atomtyping_on_gaff_paper_files()
    mol_df = load_all_gaff_paper_files()
    df = select_atomtyping()
    exit_dict = Dict{Symbol, DataFrame}
    # println((nrow(mol_df))
    for num = (1:nrow(mol_df))
        # println((num, "/", nrow(mol_df),", name: ",string(mol_df.molname[num]))
        atomtypes_list = get_atomtype(mol_df.abstract_mol[num], df)
        exit_dict = merge(exit_dict, Dict(Symbol(string(mol_df.molname[num],"atomtypes")) => atomtypes_list))
    end
    return exit_dict
end


function load_all_pdb_test_files(file_location_a::AbstractString)
    row_count = lastindex(readdir(file_location_a))
    mol_df = DataFrame([Vector{String}(undef, row_count), Vector{AbstractMolecule}(undef, row_count)], 
                        ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location_a))
        mol_df.molname[num] = string(i[1:lastindex(i)-5])
        mol_df.abstract_mol[num] = load_pdb(string(file_location_a, i))
    end
    return mol_df
end
