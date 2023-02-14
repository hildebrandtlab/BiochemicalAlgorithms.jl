@testset "Bond" begin
    # create a simple bond
    b0 = (
        idx = 0,
        a1 = 1,
        a2 = 3,
        order = BondOrder.Single,
        properties = Properties()
    )

    @test b0 isa BondTuple
    @test b0.a1 == 1
    @test b0.a2 == 3
    @test b0.order == BondOrder.Single

    # create other bonds
    for i in [3, 4, 100]
        b3 = (
            idx = 0,
            a1 = 1,
            a2 = 2,
            order = BondOrderType(i),
            properties = Properties([("bond length", i), ("is aromatic", false)])
        )

        @test b3 isa BondTuple
        @test b3.a1 == 1
        @test b3.a2 == 2
        @test !b3.properties["is aromatic"]

        if i == 3
            @test b3.order == BondOrder.Triple
            @test b3.properties["bond length"] == 3
        elseif i == 4
            @test b3.order == BondOrder.Quadruple
            @test b3.properties["bond length"] == 4
        else
            @test b3.order == BondOrder.Unknown
            @test b3.properties["bond length"] == 100
        end
    end
end
