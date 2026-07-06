function setup_optimization()
    sys = load_pdb(ball_data_path("../test/data/5PTI.pdb"))
    infer_topology!(sys)
    AmberFF(sys)
end

optimization_suite = SUITE["Optimization"]

optimization_suite["optimize_structure!"] = @benchmarkable begin
    optimize_structure!(ff; maxiters=10)
end (setup = (ff = setup_optimization()))

optimization_suite["optimize_structure_mini!"] = @benchmarkable begin
    optimize_structure_mini!(ff; epochs=100, batchsize=10)
end (setup = (ff = setup_optimization()))

optimization_suite["optimize_hydrogen_positions!"] = @benchmarkable begin
    optimize_hydrogen_positions!(ff; maxiters=10)
end (setup = (ff = setup_optimization()))