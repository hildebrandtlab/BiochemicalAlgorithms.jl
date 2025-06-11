@testitem "AminoAcidDescriptor" begin
    for T in [Float32, Float64]
        desc = AminoAcidDescriptor("AA", T(pi), T(pi), T(180))
        @test desc isa AminoAcidDescriptor{T}
        @test desc.angle_phi isa T
        @test desc.angle_psi isa T
        @test desc.angle_omega isa T
        @test desc.type isa String
        @test desc.angle_phi == T(pi)
        @test desc.angle_psi == T(pi)
        @test desc.angle_omega == T(180)
        @test desc.type == "AA"
    end

end

@testitem "PeptideBuilder" begin
    for T in [Float32, Float64]
        desc = AminoAcidDescriptor("AA", T(1), T(2), T(3))
        builder = PeptideBuilder([desc], "A", "Protein", false, FragmentDB())
        @test builder isa PeptideBuilder{T}
        @test builder.sequence isa Vector{AminoAcidDescriptor{T}}
        @test builder.chain_name isa String
        @test builder.protein_name isa String
        @test builder.is_proline isa Bool
        @test builder.fragment_db isa FragmentDB

        addAminoAcid!(builder, "A", T(1), T(2), T(3))
        @test length(builder.sequence) == 2
        @test builder.sequence[2].type == "A"

    end
end