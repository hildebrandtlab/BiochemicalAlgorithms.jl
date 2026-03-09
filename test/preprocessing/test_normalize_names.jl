@testitem "normalize_names!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

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
