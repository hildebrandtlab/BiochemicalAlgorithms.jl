@testitem "BondTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        a1 = Atom(sys, 1, Elements.H)
        a2 = Atom(sys, 2, Elements.C)
        a3 = Atom(sys, 3, Elements.O)

        b1 = Bond(a1, a2, BondOrder.Single;
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        b2 = Bond(a2, a3, BondOrder.Double)

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
        @test size(bt) == (2, 2)
        @test length(bt) == 2
        @test eltype(bt) == Bond{T}
        @test keys(bt) == [1, 2]

        # getproperty
        @test bt._sys === sys
        @test bt._idx == [b1.idx, b2.idx]

        @test bt.idx isa AbstractVector{Int}
        @test bt.idx == [b1.idx, b2.idx]
        @test bt.order isa AbstractVector{BondOrderType}
        @test bt.order == [b1.order, b2.order]

        @test bt.properties isa AbstractVector{Properties}
        @test bt.properties == [b1.properties, b2.properties]
        @test bt.flags isa AbstractVector{Flags}
        @test bt.flags == [b1.flags, b2.flags]
        @test bt.atom1_idx isa AbstractVector{Int}
        @test bt.atom1_idx == [b1.atom1_idx, b2.atom1_idx]
        @test bt.atom2_idx isa AbstractVector{Int}
        @test bt.atom2_idx == [b1.atom2_idx, b2.atom2_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(bt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :idx) == Tables.getcolumn(bt, 1) == [b1.idx, b2.idx]
        @test Tables.getcolumn(bt, :order) isa AbstractVector{BondOrderType}
        @test Tables.getcolumn(bt, 2) isa AbstractVector{BondOrderType}
        @test Tables.getcolumn(bt, :order) == Tables.getcolumn(bt, 2) == [b1.order, b2.order]

        @test Tables.getcolumn(bt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(bt, :properties) == [b1.properties, b2.properties]
        @test Tables.getcolumn(bt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(bt, :flags) == [b1.flags, b2.flags]
        @test Tables.getcolumn(bt, :atom1_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :atom1_idx) == [b1.atom1_idx, b2.atom1_idx]
        @test Tables.getcolumn(bt, :atom2_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(bt, :atom2_idx) == [b1.atom2_idx, b2.atom2_idx]

        # setproperty!
        @test_throws ErrorException bt.idx = [999, 998]
        @test_throws ErrorException bt.order = [BondOrder.Triple, BondOrder.Quadruple]

        @test_throws ErrorException bt.properties = [Properties(), Properties(:fourth => 993)]
        @test_throws ErrorException bt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException bt.atom1_idx = [997, 996]
        @test_throws ErrorException bt.atom2_idx = [995, 994]

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
        @test size(bt2) == size(bt)
        @test Tables.columnnames(bt2) == Tables.columnnames(bt)
        @test Tables.schema(bt2) == Tables.schema(bt)

        bt2 = bt[:, [:idx, :flags]]
        @test bt2 isa BondTable{T}
        @test size(bt2) == (2, 2)
        @test Tables.columnnames(bt2) == [:idx, :flags]
        @test Tables.schema(bt2).names == (:idx, :flags)
        @test Tables.schema(bt2).types == (Vector{Int}, Vector{Flags})

        bt2 = bt[2:-1:1]
        @test bt2 isa BondTable{T}
        @test length(bt2) == 2
        @test bt2[1] === bt[2]
        @test bt2[2] === bt[1]

        bt2 = bt[2:-1:1, [:idx, :flags]]
        @test bt2 isa BondTable{T}
        @test size(bt2) == (2, 2)
        @test Tables.columnnames(bt2) == [:idx, :flags]
        @test Tables.schema(bt2).names == (:idx, :flags)
        @test Tables.schema(bt2).types == (Vector{Int}, Vector{Flags})

        bt2 = bt[bt.idx .== -1]
        @test bt2 isa BondTable{T}
        @test length(bt2) == 0

        bt2 = bt[bt.idx .== b2.idx]
        @test bt2 isa BondTable{T}
        @test length(bt2) == 1
        @test only(bt2) === b2

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
        bond = Bond(atom1, atom2, BondOrder.Single)
        @test bond isa Bond{T}
        @test parent(bond) === sys
        @test parent_system(bond) === sys

        bond2 = Bond(atom2, atom3, BondOrder.Double;
            properties = Properties(:a => 1),
            flags = Flags([:A, :B])
        )
        Bond(
            Atom(sys, 1, Elements.H; frame_id = 2),
            Atom(sys, 2, Elements.C; frame_id = 2),
            BondOrder.Single
        )

        # getproperty
        @test bond.idx isa Int
        @test bond.order isa BondOrderType
        @test bond.order == BondOrder.Single

        @test bond.properties isa Properties
        @test bond.properties == Properties()
        @test bond.flags isa Flags
        @test bond.flags == Flags()
        @test bond.atom1_idx isa Int
        @test bond.atom1_idx == atom1.idx
        @test bond.atom2_idx isa Int
        @test bond.atom2_idx == atom2.idx

        @test bond._sys isa System{T}
        @test bond._idx isa Int

        @test bond2.idx isa Int
        @test bond2.order isa BondOrderType
        @test bond2.order == BondOrder.Double
        @test bond2.properties isa Properties
        @test bond2.properties == Properties(:a => 1)
        @test bond2.flags isa Flags
        @test bond2.flags == Flags([:A, :B])
        @test bond2.atom1_idx isa Int
        @test bond2.atom1_idx == atom2.idx
        @test bond2.atom2_idx isa Int
        @test bond2.atom2_idx == atom3.idx

        # setproperty!
        bond.order = BondOrder.Triple
        @test bond.order == BondOrder.Triple

        bond.properties = Properties(:first => "v1", :second => 99)
        @test length(bond.properties) == 2
        @test bond.properties[:first] == "v1"
        @test bond.properties[:second] == 99
        bond.flags = Flags([:C])
        @test length(bond.flags) == 1
        @test :C in bond.flags
        bond.atom1_idx = atom2.idx
        @test bond.atom1_idx == atom2.idx
        bond.atom2_idx = atom1.idx
        @test bond.atom2_idx == atom1.idx

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

        # nbonds
        @test nbonds(sys) isa Int
        @test nbonds(sys) == 2
        @test nbonds(sys, frame_id = 1) == 2
        @test nbonds(sys, frame_id = 2) == 1
        @test nbonds(sys, frame_id = 3) == 0
        @test nbonds(sys, frame_id = nothing) == 3

        # delete!
        @test nbonds(sys) == 2

        bidx = bond.idx
        @test delete!(bond) === nothing
        @test nbonds(sys) == 1
        @test bidx ∉ bonds(sys).idx
        @test_throws KeyError bond.idx
    end
end

@testitem "Bond/get_partner" begin
    for T in [Float32, Float64]
        sys = System{T}()
        f = Fragment(Chain(Molecule(sys)), 1)
        a1 = Atom(f, 1, Elements.C; name="C")
        a2 = Atom(f, 2, Elements.O; name="O")
        a3 = Atom(f, 3, Elements.N; name="N")
        b = Bond(a1, a2, BondOrder.Single)

        # get_partner returns the other atom of the bond
        @test get_partner(b, a1) == a2
        @test get_partner(b, a2) == a1

        # atom not in bond returns nothing
        @test isnothing(get_partner(b, a3))

        # get_partners returns both atoms
        p1, p2 = get_partners(b)
        @test p1 == a1
        @test p2 == a2
    end
end
