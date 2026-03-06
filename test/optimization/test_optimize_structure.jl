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
