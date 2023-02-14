@testset "Residue" begin
    
    r = (number = 1, 
         type = AminoAcid(1),
         chain = "A")
    @test r isa Residue
    @test r.number == 1
    @test r.type == AminoAcid(1)
    @test r.chain == "A"
end
