@testitem "MoleculeTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        m1 = Molecule(sys;
            name = "my molecule",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        m2 = Molecule(sys)

        mt = molecules(sys)

        # Tables.jl interface
        @test Tables.istable(typeof(mt))
        @test Tables.columnaccess(typeof(mt))
        @test Tables.schema(mt) isa Tables.Schema
        @test !isnothing(Tables.columns(mt))
        @test !isnothing(Tables.rows(mt))

        # AbstractArray interface
        @test size(mt) == (2, 4)
        @test length(mt) == 2
        @test eltype(mt) == Molecule{T}
        @test keys(mt) == [1, 2]

        # getproperty
        @test mt._sys === sys
        @test mt._idx == [m1.idx, m2.idx]

        @test mt.idx isa AbstractVector{Int}
        @test mt.idx == [m1.idx, m2.idx]
        @test mt.name isa AbstractVector{String}
        @test mt.name == [m1.name, m2.name]
        @test mt.properties isa AbstractVector{Properties}
        @test mt.properties == [m1.properties, m2.properties]
        @test mt.flags isa AbstractVector{Flags}
        @test mt.flags == [m1.flags, m2.flags]

        # Tables.getcolumn
        @test Tables.getcolumn(mt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, :idx) == Tables.getcolumn(mt, 1) == [m1.idx, m2.idx]
        @test Tables.getcolumn(mt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(mt, 2) isa AbstractVector{String}
        @test Tables.getcolumn(mt, :name) == Tables.getcolumn(mt, 2) == [m1.name, m2.name]
        @test Tables.getcolumn(mt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(mt, 3) isa AbstractVector{Properties}
        @test Tables.getcolumn(mt, :properties) == Tables.getcolumn(mt, 3) == [m1.properties, m2.properties]
        @test Tables.getcolumn(mt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(mt, 4) isa AbstractVector{Flags}
        @test Tables.getcolumn(mt, :flags) == Tables.getcolumn(mt, 4) == [m1.flags, m2.flags]

        # setproperty!
        @test_throws ErrorException mt.idx = [999, 998]
        @test_throws ErrorException mt.name = ["some other", "name"]
        @test_throws ErrorException mt.properties = [Properties(), Properties(:fourth => 989)]
        @test_throws ErrorException mt.flags = [Flags(), Flags([:fifth])]

        # getindex
        @test mt[1] === m1
        @test mt[2] === m2
        @test_throws BoundsError mt[0]
        @test_throws BoundsError mt[3]

        # filter
        @test filter(_ -> true, mt) == mt
        @test only(filter(m -> m.idx == m1.idx, mt)) === m1

        # collect
        mv = collect(mt)
        @test mv isa Vector{Molecule{T}}
        @test length(mv) == 2
    end
end

@testitem "Molecule" begin
    for T in [Float32, Float64]
        sys = System{T}()

        # constructors + parent
        mol = Molecule(sys)
        @test mol isa Molecule{T}
        @test parent(mol) === sys
        @test parent_system(mol) === sys

        if T == Float32
            mol_ds = Molecule()
            @test parent(mol_ds) === default_system()
            @test parent_system(mol_ds) === default_system()

            Molecule(; name = "something", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        mol2 = Molecule(sys; name = "something", properties = Properties(:a => 1), flags = Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(mol, :_row)) == 4

        # getproperty
        @test mol.idx isa Int
        @test mol.name isa String
        @test mol.name == ""
        @test mol.properties isa Properties
        @test mol.properties == Properties()
        @test mol.flags isa Flags
        @test mol.flags == Flags()

        @test mol._sys isa System{T}
        @test mol._row isa BiochemicalAlgorithms._MoleculeTableRow

        @test mol2.name == "something"
        @test mol2.properties == Properties(:a => 1)
        @test mol2.flags == Flags([:A, :B])

        # setproperty!
        mol.name = "something else"
        @test mol.name == "something else"
        mol.properties = Properties(:first => "v1", :second => 99)
        @test length(mol.properties) == 2
        @test mol.properties[:first] == "v1"
        @test mol.properties[:second] == 99
        mol.flags = Flags([:C])
        @test length(mol.flags) == 1
        @test :C in mol.flags

        # molecule_by_idx
        @test_throws KeyError molecule_by_idx(sys, -1)
        @test molecule_by_idx(sys, mol.idx) isa Molecule{T}
        @test molecule_by_idx(sys, mol.idx) == mol

        # molecules
        mv = molecules(sys)
        @test mv isa MoleculeTable{T}
        @test length(mv) == 2

        # nmolecules + push!
        @test nmolecules(sys) isa Int
        @test nmolecules(sys) == 2

        @test push!(sys, mol) === sys
        @test nmolecules(sys) == 3

        lmol = last(molecules(sys))
        @test lmol.idx != mol.idx
        @test lmol.name == mol.name
        @test lmol.properties == mol.properties
        @test lmol.flags == mol.flags

        # molecule atoms
        @test length(atoms(mol)) == 0
        @test atoms(mol) == atoms(sys, molecule_idx = mol.idx)
        @test natoms(mol) == 0
        @test natoms(mol) == natoms(sys, molecule_idx = mol.idx)

        @test Atom(mol, 1, Elements.H).molecule_idx == mol.idx
        @test length(atoms(mol)) == 1
        @test atoms(mol) == atoms(sys, molecule_idx = mol.idx)
        @test natoms(mol) == 1
        @test natoms(mol) == natoms(sys, molecule_idx = mol.idx)

        for atom in atoms(mol)
            @test parent_molecule(atom) === mol
        end
        @test parent_molecule(Atom(mol, 2, Elements.C)) === mol

        # molecule bonds
        @test length(bonds(mol)) == 0
        @test bonds(mol) == bonds(sys, molecule_idx = mol.idx)
        @test nbonds(mol) == 0
        @test nbonds(mol) == nbonds(sys, molecule_idx = mol.idx)

        Bond(mol, Atom(mol, 1, Elements.H).idx, Atom(mol, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(mol)) == 1
        @test bonds(mol) == bonds(sys, molecule_idx = mol.idx)
        @test nbonds(mol) == 1
        @test nbonds(mol) == nbonds(sys, molecule_idx = mol.idx)
    end
end
