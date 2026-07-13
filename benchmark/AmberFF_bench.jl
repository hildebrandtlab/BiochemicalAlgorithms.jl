function prepare_mol(fname)
    m = load_pdb(ball_data_path("../benchmark/data/$(fname)"))

    function atom_filter(a_idx)
        a = atom_by_idx(m, a_idx.idx)

        return !is_hetero_atom(a) && !has_flag(a, :is_deuterium)
    end

    m_filtered = copy(filter_atoms(atom_filter, m))
    infer_topology!(m_filtered)

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

# Compiled ForceField benchmarks
cff_suite = SUITE["ForceFields"]["CompiledForceField"]

cff_suite["Creation"] = @benchmarkable compile(a_ff; acc=Float32) (setup=(a_ff = AmberFF(prepare_mol("AmberFF_bench.pdb")); setup!(a_ff); update!(a_ff)))

cff = compile(a_ff; acc=Float32)

cff_suite["compute_energy!(::CompiledForceField)"] = @benchmarkable compute_energy!($cff)
cff_suite["compute_forces!(::CompiledForceField)"] = @benchmarkable compute_forces!($cff)

# Component-level benchmarks for CompiledForceField
cff_suite["_energy_stretch"] = @benchmarkable BiochemicalAlgorithms._energy_stretch($cff)
cff_suite["_energy_bend"] = @benchmarkable BiochemicalAlgorithms._energy_bend($cff)
cff_suite["_energy_torsion"] = @benchmarkable BiochemicalAlgorithms._energy_torsion($cff, $cff.proper)
cff_suite["_energy_improper"] = @benchmarkable BiochemicalAlgorithms._energy_torsion($cff, $cff.improper)

cff_suite["_sum_lj"] = @benchmarkable BiochemicalAlgorithms._sum_lj($cff.r, $cff.lj, $cff.vdw_sw, 1, length($cff.lj.i))
cff_suite["_sum_hb"] = @benchmarkable BiochemicalAlgorithms._sum_hb($cff.r, $cff.hb, $cff.vdw_sw, 1, length($cff.hb.i))
cff_suite["_sum_es"] = @benchmarkable BiochemicalAlgorithms. _sum_es($cff.r, $cff.es, $cff.es_sw, $cff.distance_dependent_dielectric, $cff.es_prefactor, 1, length($cff.es.i))


# Comparison benchmarks: reference vs compiled (Float32 acc)
comparison_suite = SUITE["ForceFields"]["Comparison"]

comparison_suite["Energy: AmberFF vs CompiledForceField (Float32)"] = @benchmarkable (compute_energy!($a_ff), compute_energy!($cff))
comparison_suite["Forces: AmberFF vs CompiledForceField (Float32)"] = @benchmarkable (compute_forces!($a_ff), compute_forces!($cff))

# Comparison benchmarks: reference vs compiled (Float64 acc)
cff64 = compile(a_ff; acc=Float64)

comparison_suite["Energy: AmberFF vs CompiledForceField (Float64)"] = @benchmarkable (compute_energy!($a_ff), compute_energy!($cff64))
comparison_suite["Forces: AmberFF vs CompiledForceField (Float64)"] = @benchmarkable (compute_forces!($a_ff), compute_forces!($cff64))
