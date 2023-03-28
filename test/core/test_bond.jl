using BiochemicalAlgorithms: _bonds

@testset "Bond" begin
    for T in [Float32, Float64]
        sys = System{T}()
        atom1 = Atom(sys, 1, Elements.H)
        atom2 = Atom(sys, 2, Elements.C)
        atom3 = Atom(sys, 3, Elements.O)

        # constructors + parent
        bond = Bond(sys, atom1.idx, atom2.idx, BondOrder.Single)
        @test bond isa Bond{T}
        @test parent(bond) === sys
        @test parent_system(bond) === sys

        if T == Float32
            bond_ds = Bond(atom1.idx, atom2.idx, BondOrder.Single)
            @test parent(bond_ds) === default_system()
            @test parent_system(bond_ds) === default_system()

            Bond(atom2.idx, atom3.idx, BondOrder.Double, Properties("a" => "b"))
        end

        bond2 = Bond(sys, atom2.idx, atom3.idx, BondOrder.Double, Properties("a" => 1))
        Bond(
            sys,
            Atom(sys, 1, Elements.H; frame_id = 2).idx,
            Atom(sys, 2, Elements.C; frame_id = 2).idx, 
            BondOrder.Single
        )

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(bond, :_row)) == 5

        # getproperty
        @test bond.idx isa Int
        @test bond.a1 isa Int
        @test bond.a1 == atom1.idx
        @test bond.a2 isa Int
        @test bond.a2 == atom2.idx
        @test bond.order isa BondOrderType
        @test bond.order == BondOrder.Single
        @test bond.properties isa Properties
        @test bond.properties == Properties()

        @test bond._sys isa System{T}
        @test bond._row isa DataFrameRow

        @test bond2.idx isa Int
        @test bond2.a1 isa Int
        @test bond2.a1 == atom2.idx
        @test bond2.a2 isa Int
        @test bond2.a2 == atom3.idx
        @test bond2.order isa BondOrderType
        @test bond2.order == BondOrder.Double
        @test bond2.properties isa Properties
        @test bond2.properties == Properties("a" => 1)

        # setproperty!
        bond.a1 = atom2.idx
        @test bond.a1 == atom2.idx
        bond.a2 = atom1.idx
        @test bond.a2 == atom1.idx
        bond.order = BondOrder.Triple
        @test bond.order == BondOrder.Triple
        bond.properties = Properties("first" => "v1", "second" => 99)
        @test length(bond.properties) == 2
        @test bond.properties["first"] == "v1"
        @test bond.properties["second"] == 99

        # bond_by_idx
        @test isnothing(bond_by_idx(sys, -1))
        @test bond_by_idx(sys, bond.idx) isa Bond{T}
        @test bond_by_idx(sys, bond.idx) == bond

        # _bonds
        df = _bonds(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(BondTuple)))
        @test copy(df[1, :]) isa BondTuple
        @test size(_bonds(sys, frame_id = 1), 1) == 2
        @test size(_bonds(sys, frame_id = 2), 1) == 1
        @test size(_bonds(sys, frame_id = 3), 1) == 0
        @test size(_bonds(sys, frame_id = nothing), 1) == 3

        # bonds_df
        df = bonds_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(BondTuple)))
        @test copy(df[1, :]) isa BondTuple
        @test size(_bonds(sys, frame_id = 1), 1) == 2
        @test size(_bonds(sys, frame_id = 2), 1) == 1
        @test size(_bonds(sys, frame_id = 3), 1) == 0
        @test size(_bonds(sys, frame_id = nothing), 1) == 3

        # bonds
        bv = bonds(sys)
        @test bv isa Vector{Bond{T}}
        @test length(bv) == 2
        @test length(bonds(sys, frame_id = 1)) == 2
        @test length(bonds(sys, frame_id = 2)) == 1
        @test length(bonds(sys, frame_id = 3)) == 0
        @test length(bonds(sys, frame_id = nothing)) == 3

        # eachbond
        @test first(eachbond(sys)) isa Bond{T}
        @test length(eachbond(sys)) == 2
        @test length(eachbond(sys, frame_id = 1)) == 2
        @test length(eachbond(sys, frame_id = 2)) == 1
        @test length(eachbond(sys, frame_id = 3)) == 0
        @test length(eachbond(sys, frame_id = nothing)) == 3

        # nbonds
        @test nbonds(sys) isa Int
        @test nbonds(sys) == 2
        @test nbonds(sys, frame_id = 1) == 2
        @test nbonds(sys, frame_id = 2) == 1
        @test nbonds(sys, frame_id = 3) == 0
        @test nbonds(sys, frame_id = nothing) == 3

        # push!
        @test push!(sys, BondTuple(
            Atom(sys, 1, Elements.H; frame_id = 3).idx,
            Atom(sys, 2, Elements.C; frame_id = 3).idx,
            BondOrder.Quadruple
        )) === sys
        @test nbonds(sys) == 2
        @test nbonds(sys, frame_id = 3) == 1
        @test nbonds(sys, frame_id = nothing) == 4
    end
end
