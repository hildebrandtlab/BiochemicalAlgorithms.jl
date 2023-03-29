@testitem "Bond Order" begin
    @test BondOrderType(1) == BondOrder.Single
    @test BondOrderType(2) == BondOrder.Double
    @test BondOrderType(3) == BondOrder.Triple
    @test BondOrderType(4) == BondOrder.Quadruple
    @test BondOrderType(50) == BondOrder.Aromatic
    @test BondOrderType(100) == BondOrder.Unknown
    @test_throws ArgumentError BondOrderType(5)
    @test_throws ArgumentError BondOrderType(-1)
end
