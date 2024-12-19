function generate_example(prepared=true)
    
    sys = load_pdb(ball_data_path("../benchmark/data/AmberFF_bench.pdb"))

    if prepared
        fdb = FragmentDB()
        normalize_names!(sys, fdb)
        reconstruct_fragments!(sys, fdb)
        build_bonds!(sys, fdb)
    end
    sys

end

kernel_iteration = SUITE["Kernel"]["Iteration"]

sys = generate_example()

kernel_iteration["AtomIteration"]=  @benchmarkable for atom in atoms($sys) end  samples = 1000000

kernel_iteration["ResidueIteration"]=  @benchmarkable for res in residues($sys) end samples = 1000000

kernel_iteration["ChainIteration"]=  @benchmarkable for chain in chains($sys) end samples = 1000000

kernel_iteration["BondIteration"]=  @benchmarkable for bond in bonds($sys) end samples = 1000000


kernel_clone = SUITE["Kernel"]["Clone"]
sys2 = generate_example(false)
kernel_clone["SystemCloning wo bonds"] = @benchmarkable s = deepcopy($sys2) samples = 1000000

sys3 = generate_example(true)
kernel_clone["SystemCloning with bonds"] = @benchmarkable s = deepcopy($sys3) samples = 1000000
