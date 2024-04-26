@testitem "AmberFF" begin
    fdb = FragmentDB()
    p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    normalize_names!(p, fdb)
    reconstruct_fragments!(p, fdb)
    build_bonds!(p, fdb)

    a_ff = AmberFF(p)

    @test compute_energy(a_ff) ≈ 1425.5979f0
    
    @test a_ff.energy["Bond Stretches"]   ≈ 1.3630637f0
    @test a_ff.energy["Angle Bends"]      ≈ 5.40766573f0
    @test a_ff.energy["Proper Torsion"]   ≈ 10.7981319f0
    @test a_ff.energy["Improper Torsion"] ≈ 3.99017335e-06
    @test a_ff.energy["Van der Waals"]    ≈ 1493.17578f0
    @test a_ff.energy["Hydrogen Bonds"]   ≈ 0f0
    @test a_ff.energy["Electrostatic"]    ≈ -85.1466827f0
end
