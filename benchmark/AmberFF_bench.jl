function prepare_mol(fname)
    fdb = FragmentDB()

    m = load_pdb(ball_data_path("../benchmark/data/$(fname)"))

    function atom_filter(a_idx)
        a = atom_by_idx(m, a_idx.idx)

        return !is_hetero_atom(a) && !has_flag(a, :is_deuterium)
    end

    m_filtered = copy(filter_atoms(atom_filter, m))

    normalize_names!(m_filtered, fdb)
    reconstruct_fragments!(m_filtered, fdb)
    build_bonds!(m_filtered, fdb)

    m_filtered
end

# benchmarking
amber_suite = SUITE["ForceFields"]["AmberFF"]

amber_suite["Creation"] = @benchmarkable AmberFF(p) (setup=(p=prepare_mol("AmberFF_bench.pdb")))

p = prepare_mol("AmberFF_bench.pdb")

amber_suite["setup!"]   = @benchmarkable setup!(a_ff) (setup=(a_ff = AmberFF($p)))

a_ff = AmberFF(p)
setup!(a_ff)

amber_suite["update!"] = @benchmarkable update!($a_ff)

update!(a_ff)

amber_suite["compute_energy!(::ForceField)"] = @benchmarkable compute_energy!($a_ff)
amber_suite["compute_forces!(::ForceField)"] = @benchmarkable compute_forces!($a_ff)

for i in a_ff.components
    amber_suite["compute_energy!(::$(i.name))"] = @benchmarkable compute_energy!($i)
end
