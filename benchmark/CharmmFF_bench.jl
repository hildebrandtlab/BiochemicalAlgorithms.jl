# benchmarking (PARAM19 with EEF1; `prepare_mol` is defined in AmberFF_bench.jl)
charmm_suite = SUITE["ForceFields"]["CharmmFF"]

# CharmmFF deletes the hydrogens of extended atoms from the given system, so every evaluation
# of the constructor needs a fresh system
charmm_suite["Creation"] = @benchmarkable CharmmFF(p) setup=(p=prepare_mol("AmberFF_bench.pdb")) evals=1

charmm_p = prepare_mol("AmberFF_bench.pdb")

# the constructor also runs setup! and update!, and removes the hydrogens from charmm_p
c_ff = CharmmFF(charmm_p)

charmm_suite["setup!"]  = @benchmarkable setup!(ff) setup=(ff = CharmmFF($charmm_p)) evals=1
charmm_suite["update!"] = @benchmarkable update!($c_ff)

charmm_suite["compute_energy!(::ForceField)"] = @benchmarkable compute_energy!($c_ff)
charmm_suite["compute_forces!(::ForceField)"] = @benchmarkable compute_forces!($c_ff)

for i in c_ff.components
    charmm_suite["compute_energy!(::$(i.name))"] = @benchmarkable compute_energy!($i)
end
