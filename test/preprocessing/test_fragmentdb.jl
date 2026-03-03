@testitem "FragmentDB" begin
    for fdb in (FragmentDB(), FragmentDB{Float32}(), FragmentDB{Float64}())
        @test length(fdb.fragments) == 33
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
            @test_broken ref.name == "Default"
            @test length(ref.atoms) == 7
            @test length(ref.bonds) == 6
        end

        # get_reference_fragment: unknown fragment returns nothing
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="UNKNOWN")
            @test isnothing(get_reference_fragment(f, fdb))
        end

        # show method
        @test contains(repr(fdb), "33 fragments")
    end
end
