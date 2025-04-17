@testitem "AminoAcidDescriptor" begin
    for T in [Float32, Float64]
        desc = AminoAcidDescriptor("AA", 1.0, 2.0, 3.0)
        @test desc isa AminoAcidDescriptor{T}
    end
end