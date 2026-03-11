@testitem "FragmentDB" begin
    using BiochemicalAlgorithms: get_reference_fragment

    for fdb in (default_fragmentdb(), FragmentDB(), FragmentDB{Float32}(), FragmentDB{Float64}())
        @test length(fdb.fragments) == 34
        @test length(fdb.name_mappings) == 6
        @test length(fdb.defaults) == 1
        @test fdb.defaults["Naming"] == "PDB"
    end

    for T in [Float32, Float64]
        # explicit path gives the same result
        fdb = FragmentDB{T}()
        fdb2 = FragmentDB{T}(ball_data_path("fragmentdb/FragmentDB.json"))
        @test fdb == fdb2

        # GLY reference fragment
        let gly = fdb.fragments["GLY"]
            @test length(gly.variants) == 4
            @test length(gly.atoms) == 10
            @test length(gly.bonds) == 9
        end

        # ALA reference fragment
        let ala = fdb.fragments["ALA"]
            @test length(ala.variants) == 4
            @test length(ala.atoms) == 13
            @test length(ala.bonds) == 12
        end

        # LYS reference fragment
        let lys = fdb.fragments["LYS"]
            @test length(lys.variants) == 4
            @test length(lys.atoms) == 25
            @test length(lys.bonds) == 24
        end

        # get_reference_fragment: known fragment returns variant
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="GLY")
            ref = get_reference_fragment(f, fdb)
            @test !isnothing(ref)
            @test ref.name == "GLY"
            @test length(ref.atoms) == 7
            @test length(ref.bonds) == 6
        end

        # get_reference_fragment: unknown fragment returns nothing
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="UNKNOWN")
            @test isnothing(get_reference_fragment(f, fdb))
        end

        # show method
        @test contains(repr(fdb), "34 fragments")
    end
end

@testitem "label_terminal_fragments!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # basic terminal labeling on AlaGlySer tripeptide
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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
            label_terminal_fragments!(sys, fdb)

            n_terms = filter(f -> :N_TERMINAL in f.flags, fragments(sys))
            c_terms = filter(f -> :C_TERMINAL in f.flags, fragments(sys))
            @test length(n_terms) == 1
            @test length(c_terms) == 1
        end
    end
end
