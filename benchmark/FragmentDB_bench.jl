function prepare_mol(fdb, steps)
    sys = load_pdb(ball_data_path("../benchmark/data/aspirin_pug_bonds.json"))
    #fdb = FragmentDB()

    if steps > 0
        normalize_names!(sys, fdb)
    end

    if steps > 1
        reconstruct_fragments!(sys, fdb)
    end

    if steps > 2
        build_bonds!(sys, fdb)
    end
    sys

end


fdb_suite = SUITE["FragmentDB"]

fdb_suite["Creation"] = @benchmarkable FragmentDB() samples = 20


fdb = FragmentDB()

fdb_suite["Normalize_Name"] = @benchmarkable normalize_names!(sys, $fdb) setup=(sys = prepare_mol($fdb,1)) samples = 50

fdb_suite["Reconstruct_Fragments"] = @benchmarkable reconstruct_fragments!(sys, $fdb) setup=(sys = prepare_mol($fdb, 2)) samples = 50

fdb_suite["Building_Bonds"] = @benchmarkable normalize_names!(sys, $fdb) (setup=(sys = prepare_mol($fdb, 3))) samples = 50
