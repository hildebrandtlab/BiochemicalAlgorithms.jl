@testitem "BondTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        a1 = Atom(sys, 1, Elements.H)
        a2 = Atom(sys, 2, Elements.C)
        a3 = Atom(sys, 3, Elements.O)

        b1 = Bond(sys, a1.idx, a2.idx, BondOrder.Single;
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        b2 = Bond(sys, a2.idx, a3.idx, BondOrder.Double)

        bt = bonds(sys)

        # AutoHashEquals, copy, and identity
        bt2 = bonds(sys)
        @test bt == bt2
        @test isequal(bt, bt2)
        @test hash(bt) == hash(bt2)
        @test bt !== bt2

        bt2 = copy(bt)
        @test bt == bt2
        @test isequal(bt, bt2)
        @test hash(bt) == hash(bt2)
        @test bt !== bt2

        # Tables.jl interface
        @test Tables.istable(typeof(bt))
        @test Tables.columnaccess(typeof(bt))
        @test Tables.schema(bt) isa Tables.Schema
        @test !isnothing(Tables.columns(bt))
        @test !isnothing(Tables.rows(bt))

        # AbstractArray interface
        @test size(bt) == (2, 4)
        @test length(bt) == 2
        @test eltype(bt) == Bond{T}
        @test keys(bt) == [1, 2]

        # getproperty
        @test bt._sys === sys
        @test bt._idx == [b1.idx, b2.idx]

        @test bt.idx isa AbstractVector{Int}
        @test bt.idx == [b1.idx, b2.idx]
        @test bt.a1 isa AbstractVector{Int}
        @test bt.a1 == [b1.a1, b2.a1]
        @test bt.a2 isa AbstractVector{Int}
        @test bt.a2 == [b1.a2, b2.a2]
        @test bt.order isa AbstractVector{BondOrderType}
        @test bt.order == [b1.order, b2.order]

        @test bt.properties isa AbstractVector{Properties}
        @test bt.properties == [b1.properties, b2.properties]
        @test bt.flags isa AbstractVector{Flags}
        @test bt.flags == [b1.flags, b2.flags]

        # Tables.getcolumn
        @test Tables.getcolumn(bt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :idx) == Tables.getcolumn(bt, 1) == [b1.idx, b2.idx]
        @test Tables.getcolumn(bt, :a1) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :a1) == Tables.getcolumn(bt, 2) == [b1.a1, b2.a1]
        @test Tables.getcolumn(bt, :a2) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, 3) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :a2) == Tables.getcolumn(bt, 3) == [b1.a2, b2.a2]
        @test Tables.getcolumn(bt, :order) isa AbstractVector{BondOrderType}
        @test Tables.getcolumn(bt, 4) isa AbstractVector{BondOrderType}
        @test Tables.getcolumn(bt, :order) == Tables.getcolumn(bt, 4) == [b1.order, b2.order]

        @test Tables.getcolumn(bt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(bt, :properties) == [b1.properties, b2.properties]
        @test Tables.getcolumn(bt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(bt, :flags) == [b1.flags, b2.flags]

        # setproperty!
        @test_throws ErrorException bt.idx = [999, 998]
        @test_throws ErrorException bt.a1 = [997, 996]
        @test_throws ErrorException bt.a2 = [995, 994]
        @test_throws ErrorException bt.order = [BondOrder.Triple, BondOrder.Quadruple]

        @test_throws ErrorException bt.properties = [Properties(), Properties(:fourth => 993)]
        @test_throws ErrorException bt.flags = [Flags(), Flags([:fifth])]

        # getindex
        @test bt[1] === b1
        @test bt[2] === b2
        @test_throws BoundsError bt[0]
        @test_throws BoundsError bt[3]

        bt2 = bt[:]
        @test bt2 isa BondTable{T}
        @test isequal(bt2, bt)
        @test bt2 == bt
        @test bt2 !== bt

        # filter
        @test filter(_ -> true, bt) == bt
        @test only(filter(b -> b.idx == b1.idx, bt)) === b1

        # collect
        bv = collect(bt)
        @test bv isa Vector{Bond{T}}
        @test length(bv) == 2

        # bonds
        @test length(bt) == 2
    end
end

@testitem "Bond" begin
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

            Bond(atom2.idx, atom3.idx, BondOrder.Double;
                properties = Properties(:a => "b"),
                flags = Flags([:A])
            )
        end

        bond2 = Bond(sys, atom2.idx, atom3.idx, BondOrder.Double;
            properties = Properties(:a => 1),
            flags = Flags([:A, :B])
        )
        Bond(
            sys,
            Atom(sys, 1, Elements.H; frame_id = 2).idx,
            Atom(sys, 2, Elements.C; frame_id = 2).idx,
            BondOrder.Single
        )

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
        @test bond.flags isa Flags
        @test bond.flags == Flags()

        @test bond._sys isa System{T}
        @test bond._idx isa Int

        @test bond2.idx isa Int
        @test bond2.a1 isa Int
        @test bond2.a1 == atom2.idx
        @test bond2.a2 isa Int
        @test bond2.a2 == atom3.idx
        @test bond2.order isa BondOrderType
        @test bond2.order == BondOrder.Double
        @test bond2.properties isa Properties
        @test bond2.properties == Properties(:a => 1)
        @test bond2.flags isa Flags
        @test bond2.flags == Flags([:A, :B])

        # setproperty!
        bond.a1 = atom2.idx
        @test bond.a1 == atom2.idx
        bond.a2 = atom1.idx
        @test bond.a2 == atom1.idx
        bond.order = BondOrder.Triple
        @test bond.order == BondOrder.Triple

        bond.properties = Properties(:first => "v1", :second => 99)
        @test length(bond.properties) == 2
        @test bond.properties[:first] == "v1"
        @test bond.properties[:second] == 99
        bond.flags = Flags([:C])
        @test length(bond.flags) == 1
        @test :C in bond.flags

        # bond_by_idx
        @test_throws KeyError bond_by_idx(sys, -1)
        @test bond_by_idx(sys, bond.idx) isa Bond{T}
        @test bond_by_idx(sys, bond.idx) == bond

        # bonds
        bv = bonds(sys)
        @test bv isa BondTable{T}
        @test length(bv) == 2
        @test length(bonds(sys, frame_id = 1)) == 2
        @test length(bonds(sys, frame_id = 2)) == 1
        @test length(bonds(sys, frame_id = 3)) == 0
        @test length(bonds(sys, frame_id = nothing)) == 3

        # nbonds + push!
        @test nbonds(sys) isa Int
        @test nbonds(sys) == 2
        @test nbonds(sys, frame_id = 1) == 2
        @test nbonds(sys, frame_id = 2) == 1
        @test nbonds(sys, frame_id = 3) == 0
        @test nbonds(sys, frame_id = nothing) == 3

        @test push!(sys, bond) === sys
        @test nbonds(sys) == 3
        @test nbonds(sys, frame_id = 1) == 3
        @test nbonds(sys, frame_id = 2) == 1
        @test nbonds(sys, frame_id = 3) == 0
        @test nbonds(sys, frame_id = nothing) == 4

        lbond = last(bonds(sys))
        @test lbond.idx != bond.idx
        @test lbond.a1 == bond.a1
        @test lbond.a2 == bond.a2
        @test lbond.order == bond.order
        @test lbond.properties == bond.properties
        @test lbond.flags == bond.flags

        # delete!
        @test nbonds(sys) == 3

        bidx = bond.idx
        @test delete!(bond) === nothing
        @test nbonds(sys) == 2
        @test bidx ∉ bonds(sys).idx
        @test_throws KeyError bond.idx
    end
end
