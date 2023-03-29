@testset "Nucleotide" begin
    # simple fragment
    n = (idx = 0, number = 1, name = "my nucleotide", properties = Properties(), flags = Flags())
    @test n isa NucleotideTuple
    @test n.number == 1
    @test n.name == "my nucleotide"
end
