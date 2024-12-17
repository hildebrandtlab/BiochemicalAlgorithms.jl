@testitem "SMARTS" begin
    for T in [Float32, Float64]
        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)

        let svec = match(SMARTSQuery(""), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test isempty(svec)
        end

        let svec = match(SMARTSQuery("C"), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test length(svec) == 8
            for at in atoms.(svec)
                @test all(e -> e == Elements.C, at.element)
            end
        end

        let svec = match(SMARTSQuery("[*H3]"), sys)
            @test svec isa Vector{<: Substructure{T}}
            @test length(svec) == 2
            @test Set(Iterators.flatten(map(sub -> atoms(sub).number, svec))) == Set([1, 7])
        end
    end
end
