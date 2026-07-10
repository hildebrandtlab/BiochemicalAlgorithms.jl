@testitem "Optimize structure" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ -374.3136f0

    # pre-minimized structure: re-optimizing should barely change energy
    e_before = compute_energy!(ff)
    @test Int(optimize_structure!(ff).retcode) == 1
    e_after = compute_energy!(ff)
    @test isapprox(e_before, e_after; atol=0.01)
end

@testitem "Optimize structure/Float64" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())

    ff = AmberFF(sys)
    @test compute_energy!(ff) isa Float64
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) < -370.0
end

@testitem "Optimize hydrogen positions" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_hydrogen_positions!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ 1421.212f0
end

@testitem "Optimize structure/minibatch" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    cff = compile(ff)
    
    e_before = compute_energy!(ff)
    @test e_before ≈ 1425.5979f0
    
    sol = optimize_structure_mini!(cff; epochs=5, batchsize=5)
    @test Int(sol.retcode) == 0
    
    e_after = compute_energy!(cff.ff)
    @test e_after < e_before
end

@testitem "Optimize structure/minibatch thread consistency" begin

    using Random

    sys1 = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys1)
    ff1 = AmberFF(sys1)
    
    Random.seed!(1234)
    cff1 = compile(ff1)
  
    
    # Run with single thread
    sys2 = deepcopy(sys1)
    ff2 = AmberFF(sys2)

    Random.seed!(1234)
    cff2 = compile(ff2, backend=:threads)

    #check both cffs should have matching components with the same seed
    @test cff1.stretch.i == cff2.stretch.i
    @test cff1.stretch.j == cff2.stretch.j
    @test cff1.bend.i == cff2.bend.i
    @test cff1.proper.i == cff2.proper.i
    @test cff1.improper.i == cff2.improper.i
    @test cff1.es.i == cff2.es.i

    r0 = collect(Float64, Iterators.flatten(atoms(cff1.ff.system).r))
    
    sol1 = optimize_structure_mini!(cff1; epochs=1, batchsize=10)
    e1 = compute_energy!(cff1.ff)
    
    r0 = collect(Float64, Iterators.flatten(atoms(cff2.ff.system).r))
    sol2 = optimize_structure_mini!(cff2; epochs=1, batchsize=10)
    e2 = compute_energy!(cff2.ff)
    @show sol1.objective, sol2.objective

    @test isapprox(sol1.u, sol2.u; rtol=1e-8)
    @test isapprox(sol1.objective, sol2.objective; rtol=1e-1)
    @test sol1.stats.iterations == sol2.stats.iterations
    @test Threads.nthreads() == 1
  
    @test isapprox(e1, e2; rtol=1e-8)  # Should be nearly identical
end

@testitem "Optimize structure/minibatch Float64" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())

    ff = AmberFF(sys)
    cff = compile(ff)
    
    sol = optimize_structure_mini!(cff; epochs=3, batchsize=8)
    @test Int(sol.retcode) == 1
    
    e_final = compute_energy!(cff.ff)
    @test e_final isa Float64
    @test e_final < 1425.0
end
