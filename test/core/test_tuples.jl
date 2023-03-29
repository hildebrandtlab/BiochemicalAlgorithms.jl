@testset "AtomTuple" begin
    # short ctor
    at = AtomTuple(1, Elements.H)
    @test at isa AtomTuple{Float32}
    @test at.idx == 0
    @test at.number == 1
    @test at.element == Elements.H
    @test at.name == ""
    @test at.atomtype == ""
    @test at.r == zeros(Vector3{Float32})
    @test at.v == zeros(Vector3{Float32})
    @test at.F == zeros(Vector3{Float32})
    @test at.formal_charge == 0
    @test at.charge == zero(Float32)
    @test at.radius == zero(Float32)
    @test !at.has_velocity
    @test !at.has_force
    @test at.properties == Properties()
    @test at.flags == Flags()

    # full ctor
    at = AtomTuple(10, Elements.C;
        idx = 9,
        name = "something",
        atomtype = "heavy atom",
        r = ones(Vector3{Float32}),
        v = 2 .* ones(Vector3{Float32}),
        F = 3 .* ones(Vector3{Float32}),
        formal_charge = 11,
        charge = Float32(12),
        radius = Float32(13),
        has_velocity = true,
        has_force = true,
        properties = Properties("a" => 14),
        flags = Flags([:A])
    )
    @test at isa AtomTuple{Float32}
    @test at.idx == 9
    @test at.number == 10
    @test at.element == Elements.C
    @test at.name == "something"
    @test at.atomtype == "heavy atom"
    @test at.r == ones(Vector3{Float32})
    @test at.v == 2 .* ones(Vector3{Float32})
    @test at.F == 3 .* ones(Vector3{Float32})
    @test at.formal_charge == 11
    @test at.charge == Float32(12)
    @test at.radius == Float32(13)
    @test at.has_velocity
    @test at.has_force
    @test at.properties == Properties("a" => 14)
    @test at.flags == Flags([:A])

    for T in [Float32, Float64]
        # short ctor
        at = AtomTuple{T}(1, Elements.H)
        @test at isa AtomTuple{T}
        @test at.idx == 0
        @test at.number == 1
        @test at.element == Elements.H
        @test at.name == ""
        @test at.atomtype == ""
        @test at.r == zeros(Vector3{T})
        @test at.v == zeros(Vector3{T})
        @test at.F == zeros(Vector3{T})
        @test at.formal_charge == 0
        @test at.charge == zero(T)
        @test at.radius == zero(T)
        @test !at.has_velocity
        @test !at.has_force
        @test at.properties == Properties()
        @test at.flags == Flags()

        # full ctor
        at = AtomTuple{T}(10, Elements.C;
            idx = 9,
            name = "something",
            atomtype = "heavy atom",
            r = ones(Vector3{T}),
            v = 2 .* ones(Vector3{T}),
            F = 3 .* ones(Vector3{T}),
            formal_charge = 11,
            charge = T(12),
            radius = T(13),
            has_velocity = true,
            has_force = true,
            properties = Properties("a" => 14),
            flags = Flags([:A])
        )
        @test at isa AtomTuple{T}
        @test at.idx == 9
        @test at.number == 10
        @test at.element == Elements.C
        @test at.name == "something"
        @test at.atomtype == "heavy atom"
        @test at.r == ones(Vector3{T})
        @test at.v == 2 .* ones(Vector3{T})
        @test at.F == 3 .* ones(Vector3{T})
        @test at.formal_charge == 11
        @test at.charge == T(12)
        @test at.radius == T(13)
        @test at.has_velocity
        @test at.has_force
        @test at.properties == Properties("a" => 14)
        @test at.flags == Flags([:A])
    end
end

@testset "BondTuple" begin
    # short ctor
    bt = BondTuple(10, 100, BondOrder.Single)
    @test bt isa BondTuple
    @test bt.idx == 0
    @test bt.a1 == 10
    @test bt.a2 == 100
    @test bt.order == BondOrder.Single
    @test bt.properties == Properties()
    @test bt.flags == Flags()

    # full ctor
    bt = BondTuple(10, 100, BondOrder.Double;
        idx = 9,
        properties = Properties("a" => 1000),
        flags = Flags([:A])
    )
    @test bt isa BondTuple
    @test bt.idx == 9
    @test bt.a1 == 10
    @test bt.a2 == 100
    @test bt.order == BondOrder.Double
    @test bt.properties == Properties("a" => 1000)
    @test bt.flags == Flags([:A])
end

@testset "MoleculeTuple" begin
    # short ctor
    mt = MoleculeTuple()
    @test mt isa MoleculeTuple
    @test mt.idx == 0
    @test mt.name == ""
    @test mt.properties == Properties()
    @test mt.flags == Flags()

    # full ctor
    mt = MoleculeTuple(;
        idx = 9,
        name = "something",
        properties = Properties("a" => 10),
        flags = Flags([:A])
    )
    @test mt isa MoleculeTuple
    @test mt.idx == 9
    @test mt.name == "something"
    @test mt.properties == Properties("a" => 10)
    @test mt.flags == Flags([:A])
end

@testset "ChainTuple" begin
    # short ctor
    ct = ChainTuple()
    @test ct isa ChainTuple
    @test ct.idx == 0
    @test ct.name == ""
    @test ct.properties == Properties()
    @test ct.flags == Flags()

    # full ctor
    ct = ChainTuple(;
        idx = 9,
        name = "something",
        properties = Properties("a" => 10),
        flags = Flags([:A])
    )
    @test ct isa ChainTuple
    @test ct.idx == 9
    @test ct.name == "something"
    @test ct.properties == Properties("a" => 10)
    @test ct.flags == Flags([:A])
end

@testset "FragmentTuple" begin
    # short ctor
    ft = FragmentTuple(1)
    @test ft isa FragmentTuple
    @test ft.idx == 0
    @test ft.number == 1
    @test ft.name == ""
    @test ft.properties == Properties()
    @test ft.flags == Flags()


    # full ctor
    ft = FragmentTuple(10;
        idx = 9,
        name = "something",
        properties = Properties("a" => 11),
        flags = Flags([:A])
    )
    @test ft isa FragmentTuple
    @test ft.idx == 9
    @test ft.number == 10
    @test ft.name == "something"
    @test ft.properties == Properties("a" => 11)
    @test ft.flags == Flags([:A])
end

@testset "NucleotideTuple" begin
    # short ctor
    nt = NucleotideTuple(1)
    @test nt isa NucleotideTuple
    @test nt.idx == 0
    @test nt.number == 1
    @test nt.name == ""
    @test nt.properties == Properties()
    @test nt.flags == Flags()

    # full ctor
    nt = NucleotideTuple(10;
        idx = 9,
        name = "something",
        properties = Properties("a" => 11),
        flags = Flags([:A])
    )
    @test nt isa NucleotideTuple
    @test nt.idx == 9
    @test nt.number == 10
    @test nt.name == "something"
    @test nt.properties == Properties("a" => 11)
    @test nt.flags == Flags([:A])
end

@testset "ResidueTuple" begin
    # short ctor
    rt = ResidueTuple(1, AminoAcid('A'))
    @test rt isa ResidueTuple
    @test rt.idx == 0
    @test rt.number == 1
    @test rt.type == AminoAcid('A')
    @test rt.properties == Properties()
    @test rt.flags == Flags()

    # full ctor
    rt = ResidueTuple(10, AminoAcid('D');
        idx = 9,
        properties = Properties("a" => 11),
        flags = Flags([:A])
    )
    @test rt isa ResidueTuple
    @test rt.idx == 9
    @test rt.number == 10
    @test rt.type == AminoAcid('D')
    @test rt.properties == Properties("a" => 11)
    @test rt.flags == Flags([:A])
end
