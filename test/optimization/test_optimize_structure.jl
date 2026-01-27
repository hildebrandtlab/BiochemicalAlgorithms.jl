@testitem "Optimize structure" tags = [:skip_ci] begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ -133.99623f0
end

@testitem "Optimize hydrogen positions" tags = [:skip_ci] begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_hydrogen_positions!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ 1424.1417f0
end
