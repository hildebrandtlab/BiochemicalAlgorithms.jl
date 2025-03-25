@testitem "FragmentDB" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        fdb = FragmentDB{T}()
        normalize_names!(sys, fdb)
        reconstruct_fragments!(sys, fdb)
        build_bonds!(sys, fdb)

        @test add_hydrogens!(sys) == 1
    end
end