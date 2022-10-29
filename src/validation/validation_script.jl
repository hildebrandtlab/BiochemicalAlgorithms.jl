using BiochemicalAlgorithms

export export_all_gaff_paper_files_to_mol2, export_all_pdb_test_files_to_mol2


function export_all_gaff_paper_files_to_mol2()
    mol_df = load_all_gaff_paper_files()
    for num = (1:nrow(mol_df))
        export_mol2(mol_df.abstract_mol[num], "../huettel-msc/export_folder/gaff_files/")
    end
end


function export_all_pdb_test_files_to_mol2(file_location_a::AbstractString, export_to_dir_location::AbstractString)
    mol_df = load_all_pdb_test_files(file_location_a)
    for num = (1:nrow(mol_df))
        export_mol2(mol_df.abstract_mol[num], export_to_dir_location)
    end
end