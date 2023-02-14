@testset "Residue" begin
    
    r = (number = 1, 
         type = AminoAcid('A'),
         chain = "A")
    @test r isa Residue
    @test r.number == 1
    @test r.type == AA_A
    @test r.chain == "A"
end
