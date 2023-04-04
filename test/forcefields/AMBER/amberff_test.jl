using BiochemicalAlgorithms

@testitem "AmberFF" begin
    fdb = FragmentDB()
    p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    normalize_names!(p, fdb)
    build_bonds!(p, fdb)
    reconstruct_fragments!(p, fdb)

    a_ff = AmberFF(p)

    @test compute_energy(a_ff) ≈ 17.5688648f0
    
    @test a_ff.energies[1] ≈ 1.3630637f0
    @test a_ff.energies[2] ≈ 5.40766573f0
    @test a_ff.energies[3] ≈ 10.798132f0
end