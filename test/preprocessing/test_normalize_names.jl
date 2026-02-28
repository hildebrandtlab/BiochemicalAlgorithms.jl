@testitem "normalize_names!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # basic normalization on AlaGlySer tripeptide
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)

            # terminal fragments should be flagged
            frags = collect(fragments(sys))
            @test :N_TERMINAL in frags[1].flags
            @test :C_TERMINAL in frags[end].flags

            # middle fragment should have neither terminal flag
            @test :N_TERMINAL ∉ frags[2].flags
            @test :C_TERMINAL ∉ frags[2].flags
        end

        # BPTI: verify terminal labeling on a real protein
        let sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
            normalize_names!(sys, fdb)

            n_terms = filter(f -> :N_TERMINAL in f.flags, fragments(sys))
            c_terms = filter(f -> :C_TERMINAL in f.flags, fragments(sys))
            @test length(n_terms) == 1
            @test length(c_terms) == 1
        end

        # system with no fragments: no crash
        let sys = System{T}()
            Molecule(sys)
            normalize_names!(sys, fdb)
            @test nfragments(sys) == 0
        end

        # invalid naming_standard: ArgumentError
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            @test_throws ArgumentError normalize_names!(sys, fdb; naming_standard="NONEXISTENT")
        end

        # verify all 6 naming schemes are present
        @test length(fdb.name_mappings) == 6
        for scheme in ["Discover", "Amber", "CHARMM", "XPLOR", "Star", "DNA"]
            @test haskey(fdb.name_mappings, scheme)
            @test fdb.name_mappings[scheme].maps_to == "PDB"
        end
    end
end
