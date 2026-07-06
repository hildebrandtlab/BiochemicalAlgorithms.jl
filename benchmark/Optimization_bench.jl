function setup_optimization(compiled_ff::Bool=false)
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)
  
    if compiled_ff
        ff = AmberFF(sys)
        return compile(ff; acc=Float32, backend=:threads)
    else
        return AmberFF(sys)
    end
end



optimization_suite = SUITE["Optimization"]

optimization_suite["optimize_structure!"] = @benchmarkable begin
    optimize_structure!(ff; maxiters=10)
end (setup = (ff = setup_optimization()))

optimization_suite["compile"] = @benchmarkable setup_optimization(true)

optimization_suite["optimize_structure_mini!"] = @benchmarkable begin
    optimize_structure_mini!(ff; epochs=100, batchsize=10)
end (setup = (ff = setup_optimization(true)))

optimization_suite["optimize_hydrogen_positions!"] = @benchmarkable begin
    optimize_hydrogen_positions!(ff; maxiters=10)
end (setup = (ff = setup_optimization()))