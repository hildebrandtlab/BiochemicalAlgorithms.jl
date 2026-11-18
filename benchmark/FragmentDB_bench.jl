function prepare_mol(fdb, steps)
    sys = load_pdb(ball_data_path("../benchmark/data/AmberFF_naming_del.pdb"))

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
fdb_suite["Creation"] = @benchmarkable FragmentDB()

fdb = FragmentDB()

fdb_suite["Normalize_Name"] = @benchmarkable normalize_names!(sys, $fdb) setup=(sys = prepare_mol($fdb,0))

fdb_suite["Reconstruct_Fragments"] = @benchmarkable reconstruct_fragments!(sys, $fdb) setup=(sys = prepare_mol($fdb, 1))

fdb_suite["Building_Bonds"] = @benchmarkable build_bonds!(sys, $fdb) setup=(sys = prepare_mol($fdb, 2))

# Larger system for reconstruct_fragments! and build_bonds!, whose cost used to
# grow quadratically with the number of atoms (#364). The test structure is
# small and already contains its hydrogens, so we remove them (to make
# reconstruct_fragments! add them back, as for a crystal structure) and append
# LARGE_COPIES copies into one system: 9 copies give ~8k atoms after
# reconstruction.
const LARGE_COPIES = 9

const large_template = let sys = load_pdb(ball_data_path("../benchmark/data/AmberFF_naming_del.pdb"))
    delete!(filter(a -> a.element == Elements.H, atoms(sys)))
    sys
end

function prepare_large_mol(fdb, template, steps)
    sys = deepcopy(template)
    append!(sys, fill(template, LARGE_COPIES - 1)...)

    normalize_names!(sys, fdb)

    if steps > 1
        reconstruct_fragments!(sys, fdb)
    end
    sys
end

# These take seconds per run before the #364 fixes, so cap the number of samples.
fdb_suite["Reconstruct_Fragments_Large"] = @benchmarkable reconstruct_fragments!(sys, $fdb) setup=(sys = prepare_large_mol($fdb, $large_template, 1)) evals=1 samples=5

fdb_suite["Building_Bonds_Large"] = @benchmarkable build_bonds!(sys, $fdb) setup=(sys = prepare_large_mol($fdb, $large_template, 2)) evals=1 samples=5
