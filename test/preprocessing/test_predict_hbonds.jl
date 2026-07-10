@testitem "predict_hbonds!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # KABSCH_SANDER on 2ptc (h_bond_from_donor=true, default)
        let sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"), T)
            infer_topology!(sys, fdb)

            predict_hbonds!(sys, :KABSCH_SANDER)
            hbonds = hydrogen_bonds(sys)
            @test length(hbonds) == 181

            # all Kabsch-Sander H-bonds are backbone bonds
            bb = backbone_hydrogen_bonds(sys)
            @test length(bb) == 181

            # verify H-bond flags
            for b in hbonds
                @test :TYPE__HYDROGEN in b.flags
            end
        end

        # KABSCH_SANDER with h_bond_from_donor=false
        let sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"), T)
            infer_topology!(sys, fdb)

            predict_hbonds!(sys, :KABSCH_SANDER, false)
            @test length(hydrogen_bonds(sys)) == 175
        end

        # KABSCH_SANDER on bpti
        let sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
            infer_topology!(sys, fdb)

            predict_hbonds!(sys, :KABSCH_SANDER)
            @test length(hydrogen_bonds(sys)) == 31
            @test length(backbone_hydrogen_bonds(sys)) == 31
        end

        # unknown method: ErrorException
        let sys = System{T}()
            @test_throws ErrorException predict_hbonds!(sys, :NONEXISTENT)
        end

        # idempotency: predicting twice should not duplicate bonds
        let sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
            infer_topology!(sys, fdb)

            predict_hbonds!(sys, :KABSCH_SANDER)
            n = length(hydrogen_bonds(sys))
            @test n > 0

            predict_hbonds!(sys, :KABSCH_SANDER)
            @test length(hydrogen_bonds(sys)) == n
        end
    end
end
