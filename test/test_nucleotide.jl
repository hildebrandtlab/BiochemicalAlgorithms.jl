@testset "Nucleotide" begin
    # simple fragment
    n = ( number = 1, name = "my nucleotide", chain = "chain")
    @test n isa Nucleotide
    @test n.number == 1
    @test n.name == "my nucleotide"
    @test n.chain == "chain"

    # incomplete fragment
    n2 = ( number = 1, name = "my nucleotide")
    @test !isa(n2, Nucleotide)

    # nucleotide is compatible to fragment
    @test_broken !isa(n, Fragment)
end
