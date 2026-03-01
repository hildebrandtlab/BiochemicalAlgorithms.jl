@testitem "SMARTS" begin
    for T in [Float32, Float64]
        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)

        # empty pattern → no matches
        let svec = match(SMARTSQuery(""), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test isempty(svec)
        end

        # carbon atoms
        let svec = match(SMARTSQuery("C"), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test length(svec) == 8
            for at in atoms.(svec)
                @test all(e -> e == Elements.C, at.element)
            end
        end

        # hydrogen count query
        let svec = match(SMARTSQuery("[*H3]"), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test length(svec) == 2
            @test Set(Iterators.flatten(map(sub -> atoms(sub).number, svec))) == Set([1, 7])
        end

        # nitrogen atoms
        let svec = match(SMARTSQuery("N"), sys)
            @test length(svec) == 3
            for at in atoms.(svec)
                @test all(e -> e == Elements.N, at.element)
            end
        end

        # oxygen atoms
        let svec = match(SMARTSQuery("O"), sys)
            @test length(svec) == 5
            for at in atoms.(svec)
                @test all(e -> e == Elements.O, at.element)
            end
        end

        # atomic number query (carbon = 6)
        let svec = match(SMARTSQuery("[#6]"), sys)
            @test length(svec) == 8
        end

        # bond patterns: C-C, C-O, C-N
        @test length(match(SMARTSQuery("CC"), sys)) == 10
        @test length(match(SMARTSQuery("CO"), sys)) == 3
        @test length(match(SMARTSQuery("CN"), sys)) == 5

        # N with one hydrogen
        @test length(match(SMARTSQuery("[NH]"), sys)) == 2

        # O with one hydrogen
        @test length(match(SMARTSQuery("[OH]"), sys)) == 3
    end
end
