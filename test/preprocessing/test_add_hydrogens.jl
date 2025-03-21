@testitem "FragmentDB" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    @test add_hydrogens!(sys) == 1
end