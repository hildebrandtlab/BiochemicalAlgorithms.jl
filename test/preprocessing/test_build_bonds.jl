@testitem "build_bonds!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            empty!(sys._bonds)

            normalize_names!(sys, fdb)
            build_bonds!(sys, fdb)
            @test nbonds.(fragments(sys)) == [12, 8, 12]
        end

        let sys = load_pdb(ball_data_path("../test/data/AmberFF_bench.pdb"), T)
            empty!(sys._bonds)

            normalize_names!(sys, fdb)
            build_bonds!(sys, fdb)
            @test nbonds(sys) == 906
        end

        let sys = load_pdb(ball_data_path("../test/data/ACE_test_A.pdb"), T)
            normalize_names!(sys, fdb)
            build_bonds!(sys, fdb)
            @test nbonds(sys) == 1666
        end

        let sys = load_pdb(ball_data_path("../test/data/ACE_test_B.pdb"), T)
            normalize_names!(sys, fdb)
            build_bonds!(sys, fdb)
            @test nbonds(sys) == 468
        end
    end
end
