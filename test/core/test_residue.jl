@testitem "Residue" begin
    using BioSymbols
    
    r = (
        idx = 0,
        number = 1, 
        type = AminoAcid('A'),
        properties = Properties()
    )
    @test r isa ResidueTuple
    @test r.number == 1
    @test r.type == AA_A
end
