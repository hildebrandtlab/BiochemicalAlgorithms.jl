@testitem "Atom" begin
    using DataFrames

    using BiochemicalAlgorithms: _atoms, _bonds

    for T in [Float32, Float64]
        at = AtomTuple{T}(1, Elements.H;
            name = "my fancy atom",
            atom_type = "heavy",
            r = Vector3{T}(1, 2, 4),
            v = Vector3{T}(1, 1, 1),
            formal_charge = 1,
            charge = T(0.2),
            radius = T(1.01)
        )
        sys = System{T}()

        # constructors + parent
        atom = Atom(sys, at)
        @test atom isa Atom{T}
        @test parent(atom) === sys
        @test parent_system(atom) === sys
        T == Float32 && @test Atom(at.number, at.element) isa Atom{T}
        T == Float32 && @test parent(Atom(at)) === default_system()
        T == Float32 && @test parent_system(Atom(at)) === default_system()

        atom2 = Atom(sys, at, frame_id = 10, molecule_id = 11, chain_id = 12, fragment_id = 13, 
            nucleotide_id = 14, residue_id = 15)

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(atom, :_row)) == 13

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
        @test atom.r isa Position{T}
        @test atom.r == at.r
        @test atom.v isa Velocity{T}
        @test atom.v == at.v
        @test atom.F isa Force{T}
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

        @test atom._sys isa System{T}
        @test atom._row isa AtomTableRow{T}

        @test_throws ErrorException atom.frame_id
        @test_throws ErrorException atom.molecule_id
        @test_throws ErrorException atom.chain_id
        @test_throws ErrorException atom.fragment_id
        @test_throws ErrorException atom.nucleotide_id
        @test_throws ErrorException atom.residue_id

        @test atom._row.frame_id isa Int
        @test atom._row.frame_id == 1
        @test isnothing(atom._row.molecule_id)
        @test isnothing(atom._row.chain_id)
        @test isnothing(atom._row.fragment_id)
        @test isnothing(atom._row.nucleotide_id)
        @test isnothing(atom._row.residue_id)

        @test atom2._row.frame_id isa Int
        @test atom2._row.frame_id == 10
        @test atom2._row.molecule_id isa Int
        @test atom2._row.molecule_id == 11
        @test atom2._row.chain_id isa Int
        @test atom2._row.chain_id == 12
        @test atom2._row.fragment_id isa Int
        @test atom2._row.fragment_id == 13
        @test atom2._row.nucleotide_id isa Int
        @test atom2._row.nucleotide_id == 14
        @test atom2._row.residue_id isa Int
        @test atom2._row.residue_id == 15

        # setproperty!
        atom.number = 42
        @test atom.number == 42
        atom.element = Elements.C
        @test atom.element == Elements.C
        atom.name = "another name"
        @test atom.name == "another name"
        atom.atom_type = "none"
        @test atom.atom_type == "none"
        atom.r = Position(T[10, 20, 30])
        @test atom.r == Position(T[10, 20, 30])
        atom.v = Velocity(T[100, 200, 300])
        @test atom.v == Velocity(T[100, 200, 300])
        atom.F = Force(T[1000, 2000, 3000])
        @test atom.F == Force(T[1000, 2000, 3000])
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

        @test_throws ErrorException atom.frame_id = 0
        @test_throws ErrorException atom.molecule_id = 0
        @test_throws ErrorException atom.chain_id = 0
        @test_throws ErrorException atom.fragment_id = 0
        @test_throws ErrorException atom.nucleotide_id = 0
        @test_throws ErrorException atom.residue_id = 0

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
        atom3 = Atom(mol, at.number, at.element, at.name; frame_id = 10)
        @test atom_by_name(mol, atom3.name; frame_id = 10) == atom3

        # atoms
        avec = atoms(sys)
        @test avec isa Vector{Atom{T}}
        @test length(avec) == 1
        @test length(atoms(sys, frame_id = 1)) == 1
        @test length(atoms(sys, frame_id = 2)) == 0
        @test length(atoms(sys, frame_id = 10)) == 2
        @test length(atoms(sys, frame_id = nothing)) == 3
        @test length(atoms(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15)) == 1

        # atoms_df
        df = atoms_df(sys)
        @test df isa DataFrame
        @test size(df) == (1, length(fieldnames(AtomTuple{T})))
        @test copy(df[1, :]) isa AtomTuple{T}
        @test size(atoms_df(sys, frame_id = 1), 1) == 1
        @test size(atoms_df(sys, frame_id = 2), 1) == 0
        @test size(atoms_df(sys, frame_id = 10), 1) == 2
        @test size(atoms_df(sys, frame_id = nothing), 1) == 3
        @test size(atoms_df(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15), 1) == 1

        # eachatom
        @test length(collect(eachatom(sys))) == 1
        @test first(eachatom(sys)) isa Atom{T}
        @test length(collect(eachatom(sys, frame_id = 1))) == 1
        @test length(collect(eachatom(sys, frame_id = 2))) == 0
        @test length(collect(eachatom(sys, frame_id = 10))) == 2
        @test length(collect(eachatom(sys, frame_id = nothing))) == 3
        @test length(collect(eachatom(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15))) == 1

        # natoms + push!
        @test natoms(sys) isa Int
        @test natoms(sys) == 1
        @test natoms(sys, frame_id = 1) == 1
        @test natoms(sys, frame_id = 2) == 0
        @test natoms(sys, frame_id = 10) == 2
        @test natoms(sys, frame_id = nothing) == 3
        @test natoms(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15) == 1

        @test push!(sys, at) === sys
        @test natoms(sys) == 2
        @test push!(sys, at, frame_id = 100, molecule_id = 101, chain_id = 102, fragment_id = 103, 
            nucleotide_id = 104, residue_id = 105) === sys
        @test natoms(sys) == 2
        @test natoms(sys, frame_id = 100) == 1

        # test is_geminal, is_vicinal
        a = Atom(sys, at)
        b = Atom(sys, at)
        c = Atom(sys, at)

        @test !is_geminal(a, a)
        @test !is_geminal(a, b)

        Bond(sys, a.idx, b.idx, BondOrder.Single)
        Bond(sys, b.idx, c.idx, BondOrder.Single)
        
        @test !is_geminal(a, b)
        @test is_geminal(a, c)
        @test is_geminal(c, a)

        a = Atom(sys, at)
        b = Atom(sys, at)
        c = Atom(sys, at)
        d = Atom(sys, at)

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
        @test size(_bonds(atom), 1) == 0
        @test size(bonds_df(atom), 1) == 0
        @test length(bonds(atom)) == 0
        @test length(eachbond(atom)) == 0
        @test nbonds(atom) == 0

        @test push!(sys, BondTuple(
            atom.idx,
            Atom(sys, 2, Elements.C).idx,
            BondOrder.Single
        )) === sys
        @test size(_bonds(atom), 1) == 1
        @test size(bonds_df(atom), 1) == 1
        @test length(bonds(atom)) == 1
        @test length(eachbond(atom)) == 1
        @test nbonds(atom) == 1

        @test push!(sys, BondTuple(
            Atom(sys, 3, Elements.C).idx,
            atom.idx,
            BondOrder.Double
        )) === sys
        @test size(_bonds(atom), 1) == 2
        @test size(bonds_df(atom), 1) == 2
        @test length(bonds(atom)) == 2
        @test length(eachbond(atom)) == 2
        @test nbonds(atom) == 2

        @test push!(sys, BondTuple(
            Atom(sys, 4, Elements.C).idx,
            Atom(sys, 5, Elements.C).idx,
            BondOrder.Double
        )) === sys
        @test size(_bonds(atom), 1) == 2
        @test size(bonds_df(atom), 1) == 2
        @test length(bonds(atom)) == 2
        @test length(eachbond(atom)) == 2
        @test nbonds(atom) == 2
    end
end
