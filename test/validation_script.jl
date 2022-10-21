using BiochemicalAlgorithms: load_all, export_mol2

export export_all_to_mol2


function export_all_to_mol2()
    mol_df = load_all()
    for num = (1:nrow(mol_df))
        export_mol2(mol_df.abstract_mol[num])
    end
end
