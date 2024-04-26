@testitem "AtomTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        a1 = Atom(sys, 1, Elements.H;
            name = "my atom",
            atom_type = "my atom type",
            r = ones(Vector3{T}),
            v = 2 .* ones(Vector3{T}),
            F = 3 .* ones(Vector3{T}),
            formal_charge = 4,
            charge = T(5),
            radius = T(6),
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        a2 = Atom(sys, 2, Elements.C)

        at = atoms(sys)
        
        # Tables.jl interface
        @test Tables.istable(typeof(at))
        @test Tables.columnaccess(typeof(at))
        @test Tables.schema(at) isa Tables.Schema
        @test !isnothing(Tables.columns(at))
        @test !isnothing(Tables.rows(at))

        # AbstractArray interface
        @test size(at) == (2, 11)
        @test length(at) == 2
        @test eltype(at) == Atom{T}
        @test keys(at) == [1, 2]

        # getproperty
        @test at._sys === sys
        @test at._idx == [a1.idx, a2.idx]

        @test at.idx isa AbstractVector{Int}
        @test at.idx == [a1.idx, a2.idx]
        @test at.number isa AbstractVector{Int}
        @test at.number == [a1.number, a2.number]
        @test at.element isa AbstractVector{ElementType}
        @test at.element == [a1.element, a2.element]
        @test at.name isa AbstractVector{String}
        @test at.name == [a1.name, a2.name]
        @test at.atom_type isa AbstractVector{String}
        @test at.atom_type == [a1.atom_type, a2.atom_type]
        @test at.r isa AbstractVector{Vector3{T}}
        @test at.r == [a1.r, a2.r]
        @test at.v isa AbstractVector{Vector3{T}}
        @test at.v == [a1.v, a2.v]
        @test at.F isa AbstractVector{Vector3{T}}
        @test at.F == [a1.F, a2.F]
        @test at.formal_charge isa AbstractVector{Int}
        @test at.formal_charge == [a1.formal_charge, a2.formal_charge]
        @test at.charge isa AbstractVector{T}
        @test at.charge == [a1.charge, a2.charge]
        @test at.radius isa AbstractVector{T}
        @test at.radius == [a1.radius, a2.radius]

        @test at.properties isa AbstractVector{Properties}
        @test at.properties == [a1.properties, a2.properties]
        @test at.flags isa AbstractVector{Flags}
        @test at.flags == [a1.flags, a2.flags]
        @test at.frame_id isa AbstractVector{Int}
        @test at.frame_id == [a1.frame_id, a2.frame_id]
        @test at.molecule_idx isa AbstractVector{MaybeInt}
        @test at.molecule_idx == [a1.molecule_idx, a2.molecule_idx]
        @test at.chain_idx isa AbstractVector{MaybeInt}
        @test at.chain_idx == [a1.chain_idx, a2.chain_idx]
        @test at.fragment_idx isa AbstractVector{MaybeInt}
        @test at.fragment_idx == [a1.fragment_idx, a2.fragment_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(at, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(at, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(at, :idx) == Tables.getcolumn(at, 1) == [a1.idx, a2.idx]
        @test Tables.getcolumn(at, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(at, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(at, :number) == Tables.getcolumn(at, 2) == [a1.number, a2.number]
        @test Tables.getcolumn(at, :element) isa AbstractVector{ElementType}
        @test Tables.getcolumn(at, 3) isa AbstractVector{ElementType}
        @test Tables.getcolumn(at, :element) == Tables.getcolumn(at, 3) == [a1.element, a2.element]
        @test Tables.getcolumn(at, :name) isa AbstractVector{String}
        @test Tables.getcolumn(at, 4) isa AbstractVector{String}
        @test Tables.getcolumn(at, :name) == Tables.getcolumn(at, 4) == [a1.name, a2.name]
        @test Tables.getcolumn(at, :atom_type) isa AbstractVector{String}
        @test Tables.getcolumn(at, 5) isa AbstractVector{String}
        @test Tables.getcolumn(at, :atom_type) == Tables.getcolumn(at, 5) == [a1.atom_type, a2.atom_type]
        @test Tables.getcolumn(at, :r) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, 6) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, :r) == Tables.getcolumn(at, 6) == [a1.r, a2.r]
        @test Tables.getcolumn(at, :v) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, 7) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, :v) == Tables.getcolumn(at, 7) == [a1.v, a2.v]
        @test Tables.getcolumn(at, :F) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, 8) isa AbstractVector{Vector3{T}}
        @test Tables.getcolumn(at, :F) == Tables.getcolumn(at, 8) == [a1.F, a2.F]
        @test Tables.getcolumn(at, :formal_charge) isa AbstractVector{Int}
        @test Tables.getcolumn(at, 9) isa AbstractVector{Int}
        @test Tables.getcolumn(at, :formal_charge) == Tables.getcolumn(at, 9) == [a1.formal_charge, a2.formal_charge]
        @test Tables.getcolumn(at, :charge) isa AbstractVector{T}
        @test Tables.getcolumn(at, 10) isa AbstractVector{T}
        @test Tables.getcolumn(at, :charge) == Tables.getcolumn(at, 10) == [a1.charge, a2.charge]
        @test Tables.getcolumn(at, :radius) isa AbstractVector{T}
        @test Tables.getcolumn(at, 11) isa AbstractVector{T}
        @test Tables.getcolumn(at, :radius) == Tables.getcolumn(at, 11) == [a1.radius, a2.radius]

        @test Tables.getcolumn(at, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(at, :properties) == [a1.properties, a2.properties]
        @test Tables.getcolumn(at, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(at, :flags) == [a1.flags, a2.flags]
        @test Tables.getcolumn(at, :frame_id) isa AbstractVector{Int}
        @test Tables.getcolumn(at, :frame_id) == [a1.frame_id, a2.frame_id]
        @test Tables.getcolumn(at, :molecule_idx) isa AbstractVector{MaybeInt}
        @test Tables.getcolumn(at, :molecule_idx) == [a1.molecule_idx, a2.molecule_idx]
        @test Tables.getcolumn(at, :chain_idx) isa AbstractVector{MaybeInt}
        @test Tables.getcolumn(at, :chain_idx) == [a1.chain_idx, a2.chain_idx]
        @test Tables.getcolumn(at, :fragment_idx) isa AbstractVector{MaybeInt}
        @test Tables.getcolumn(at, :fragment_idx) == [a1.fragment_idx, a2.fragment_idx]

        # setproperty!
        @test_throws ErrorException at.idx = [999, 998]
        @test_throws ErrorException at.number = [997, 996]
        @test_throws ErrorException at.element = [Elements.S, Elements.N]
        @test_throws ErrorException at.name = ["some other", "names"]
        @test_throws ErrorException at.atom_type = ["some other", "types"]
        @test_throws ErrorException at.r = [zeros(Vector3{T}), ones(Vector3{T})]
        @test_throws ErrorException at.v = [zeros(Vector3{T}), ones(Vector3{T})]
        @test_throws ErrorException at.F = [zeros(Vector3{T}), ones(Vector3{T})]
        @test_throws ErrorException at.formal_charge = [995, 994]
        @test_throws ErrorException at.charge = T[993, 992]
        @test_throws ErrorException at.radius = T[991, 990]

        @test_throws ErrorException at.properties = [Properties(), Properties(:fourth => 989)]
        @test_throws ErrorException at.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException at.frame_id = [988, 987]
        @test_throws ErrorException at.molecule_idx = [986, 985]
        @test_throws ErrorException at.chain_idx = [984, 983]
        @test_throws ErrorException at.fragment_idx = [982, 981]

        # getindex
        @test at[1] === a1
        @test at[2] === a2
        @test_throws BoundsError at[0]
        @test_throws BoundsError at[3]

        # filter
        @test filter(_ -> true, at) == at
        @test only(filter(a -> a.idx == a1.idx, at)) === a1

        # collect
        av = collect(at)
        @test av isa Vector{Atom{T}}
        @test length(av) == 2
    end
end

@testitem "Atom" begin
    for T in [Float32, Float64]
        at = (
            number = 1,
            element = Elements.H,
            name = "my atom",
            atom_type = "my atom type",
            r = ones(Vector3{T}),
            v = 2 .* ones(Vector3{T}),
            F = 3 .* ones(Vector3{T}),
            formal_charge = 4,
            charge = T(5),
            radius = T(6),
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        sys = System{T}()

        # constructors + parent
        atom = Atom(sys, at.number, at.element;
            name = at.name,
            atom_type = at.atom_type,
            r = at.r,
            v = at.v,
            F = at.F,
            formal_charge = at.formal_charge,
            charge = at.charge,
            radius = at.radius,
            properties = at.properties,
            flags = at.flags
        )
        @test atom isa Atom{T}
        @test parent(atom) === sys
        @test parent_system(atom) === sys
        T == Float32 && @test parent(Atom(at.number, at.element)) === default_system()
        T == Float32 && @test parent_system(Atom(at.number, at.element)) === default_system()

        @test isnothing(parent_molecule(atom))
        @test isnothing(parent_chain(atom))
        @test isnothing(parent_fragment(atom))

        atom2 = Atom(sys, at.number, at.element;
            name = at.name,
            atom_type = at.atom_type,
            r = at.r,
            v = at.v,
            formal_charge = at.formal_charge,
            charge = at.charge,
            radius = at.radius,
            frame_id = 10,
            molecule_idx = 11,
            chain_idx = 12,
            fragment_idx = 13
        )

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(atom._row) == 11

        # getproperty
        @test atom.idx isa Int
        @test atom.number isa Int
        @test atom.number == at.number
        @test atom.element isa ElementType
        @test atom.element == at.element
        @test atom.name isa String
        @test atom.name == at.name
        @test atom.atom_type isa String
        @test atom.atom_type == at.atom_type
        @test atom.r isa Vector3{T}
        @test atom.r == at.r
        @test atom.v isa Vector3{T}
        @test atom.v == at.v
        @test atom.F isa Vector3{T}
        @test atom.F == at.F
        @test atom.formal_charge isa Int
        @test atom.formal_charge == at.formal_charge
        @test atom.charge isa T
        @test atom.charge == at.charge
        @test atom.radius isa T
        @test atom.radius == at.radius

        @test atom.properties isa Properties
        @test atom.properties == at.properties
        @test atom.flags isa Flags
        @test atom.flags == at.flags
        @test atom.frame_id isa Int
        @test atom.frame_id == 1
        @test isnothing(atom.molecule_idx)
        @test isnothing(atom.chain_idx)
        @test isnothing(atom.fragment_idx)

        @test atom._sys isa System{T}
        @test atom._row isa BiochemicalAlgorithms._AtomTableRow{T}

        @test atom2.frame_id isa Int
        @test atom2.frame_id == 10
        @test atom2.molecule_idx isa Int
        @test atom2.molecule_idx == 11
        @test atom2.chain_idx isa Int
        @test atom2.chain_idx == 12
        @test atom2.fragment_idx isa Int
        @test atom2.fragment_idx == 13

        # setproperty!
        atom.number = 42
        @test atom.number == 42
        atom.element = Elements.C
        @test atom.element == Elements.C
        atom.name = "another name"
        @test atom.name == "another name"
        atom.atom_type = "none"
        @test atom.atom_type == "none"
        atom.r = Vector3{T}(10, 20, 30)
        @test atom.r == Vector3{T}(10, 20, 30)
        atom.v = Vector3{T}(100, 200, 300)
        @test atom.v == Vector3{T}(100, 200, 300)
        atom.F = Vector3{T}(1000, 2000, 3000)
        @test atom.F == Vector3{T}(1000, 2000, 3000)
        atom.formal_charge = 2
        @test atom.formal_charge == 2
        atom.charge = -one(T)
        @test atom.charge == -one(T)
        atom.radius = one(T) / 2
        @test atom.radius == one(T) / 2

        atom.properties = Properties(:first => "v1", :second => 99)
        @test length(atom.properties) == 2
        @test atom.properties[:first] == "v1"
        @test atom.properties[:second] == 99
        atom.flags = Flags([:A, :B])
        @test length(atom.flags) == 2
        @test :A in atom.flags
        @test :B in atom.flags

        atom3 = Atom(System{T}(), at.number, at.element; name = at.name, frame_id = 10)
        atom3.frame_id = 999
        @test atom3.frame_id == 999
        atom3.molecule_idx = 998
        @test atom3.molecule_idx == 998
        atom3.chain_idx = 997
        @test atom3.chain_idx == 997
        atom3.fragment_idx = 996
        @test atom3.fragment_idx == 996

        # atom_by_idx
        @test_throws KeyError atom_by_idx(sys, -1)
        @test atom_by_idx(sys, atom.idx) isa Atom{T}
        @test atom_by_idx(sys, atom.idx) == atom

        # atom_by_name
        @test isnothing(atom_by_name(sys, "invalid"))
        @test atom_by_name(sys, atom.name) isa Atom{T}
        @test atom_by_name(sys, atom.name) == atom
        @test atom_by_name(sys, atom.name; frame_id = 1) == atom
        @test isnothing(atom_by_name(sys, atom.name; frame_id = 9999))
        @test isnothing(atom_by_name(sys, atom2.name))
        @test atom_by_name(sys, atom2.name; frame_id = 10) == atom2
        mol = Molecule(sys)
        atom3 = Atom(mol, at.number, at.element; name = at.name, frame_id = 10)
        @test atom_by_name(mol, atom3.name; frame_id = 10) == atom3

        # atoms
        avec = atoms(sys)
        @test avec isa AtomTable{T}
        @test length(avec) == 1
        @test length(atoms(sys, frame_id = 1)) == 1
        @test length(atoms(sys, frame_id = 2)) == 0
        @test length(atoms(sys, frame_id = 10)) == 2
        @test length(atoms(sys, frame_id = nothing)) == 3
        @test length(atoms(sys, frame_id = nothing, molecule_idx =11, chain_idx = 12, fragment_idx = 13)) == 1

        # natoms + push!
        @test natoms(sys) isa Int
        @test natoms(sys) == 1
        @test natoms(sys, frame_id = 1) == 1
        @test natoms(sys, frame_id = 2) == 0
        @test natoms(sys, frame_id = 10) == 2
        @test natoms(sys, frame_id = nothing) == 3
        @test natoms(sys, frame_id = nothing, molecule_idx =11, chain_idx = 12, fragment_idx = 13) == 1

        @test push!(sys, atom) === sys
        @test natoms(sys) == 2
        @test push!(sys, atom, frame_id = 100, molecule_idx = 101, chain_idx = 102, fragment_idx = 103) === sys
        @test natoms(sys) == 2
        @test natoms(sys, frame_id = 100) == 1

        latom = last(atoms(sys))
        @test latom.idx != atom.idx
        @test latom.number == atom.number
        @test latom.element == atom.element
        @test latom.name == atom.name
        @test latom.atom_type == atom.atom_type
        @test latom.r == atom.r
        @test latom.v == atom.v
        @test latom.F == atom.F
        @test latom.formal_charge == atom.formal_charge
        @test latom.charge == atom.charge
        @test latom.radius == atom.radius
        @test latom.properties == atom.properties
        @test latom.flags == atom.flags
        @test latom.frame_id == atom.frame_id
        @test latom.molecule_idx == atom.molecule_idx
        @test latom.chain_idx == atom.chain_idx
        @test latom.fragment_idx == atom.fragment_idx

        # test is_geminal, is_vicinal
        a = Atom(sys, at.number, at.element)
        b = Atom(sys, at.number, at.element)
        c = Atom(sys, at.number, at.element)

        @test !is_geminal(a, a)
        @test !is_geminal(a, b)

        Bond(sys, a.idx, b.idx, BondOrder.Single)
        Bond(sys, b.idx, c.idx, BondOrder.Single)
        
        @test !is_geminal(a, b)
        @test is_geminal(a, c)
        @test is_geminal(c, a)

        a = Atom(sys, at.number, at.element)
        b = Atom(sys, at.number, at.element)
        c = Atom(sys, at.number, at.element)
        d = Atom(sys, at.number, at.element)

        @test !is_vicinal(a, a)
        @test !is_vicinal(a, b)

        Bond(sys, a.idx, b.idx, BondOrder.Single)
        Bond(sys, b.idx, c.idx, BondOrder.Single)
        Bond(sys, c.idx, d.idx, BondOrder.Single)

        @test !is_vicinal(a, a)
        @test !is_vicinal(a, c)
        @test is_vicinal(a, d)
        @test is_vicinal(d, a)

        # atom bonds
        @test length(bonds(atom)) == 0
        @test nbonds(atom) == 0

        @test parent(Bond(
            sys,
            atom.idx,
            Atom(sys, 2, Elements.C).idx,
            BondOrder.Single
        )) === sys
        @test length(bonds(atom)) == 1
        @test nbonds(atom) == 1

        @test parent(Bond(
            sys,
            Atom(sys, 3, Elements.C).idx,
            atom.idx,
            BondOrder.Double
        )) === sys
        @test length(bonds(atom)) == 2
        @test nbonds(atom) == 2

        @test parent(Bond(
            sys,
            Atom(sys, 4, Elements.C).idx,
            Atom(sys, 5, Elements.C).idx,
            BondOrder.Double
        )) === sys
        @test length(bonds(atom)) == 2
        @test nbonds(atom) == 2
    end
end
