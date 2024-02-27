@testitem "Molecule" begin
    using BiochemicalAlgorithms: _AtomTable, _BondTable, _MoleculeTableRow, _molecules, _atoms, _bonds
    using DataFrames

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

            Molecule("something", Properties(:a => "b"), Flags([:A]))
        end

        mol2 = Molecule(sys, "something", Properties(:a => 1), Flags([:A, :B]))

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
        @test mol._row isa _MoleculeTableRow

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

        # eachmolecule
        @test first(eachmolecule(sys)) isa Molecule{T}
        @test length(eachmolecule(sys)) == 2

        # nmolecules
        @test nmolecules(sys) isa Int
        @test nmolecules(sys) == 2

        # molecule atoms
        @test length(collect(_atoms(mol))) == 0
        @test Tables.materializer(_AtomTable{T})(_atoms(mol)) == 
            Tables.materializer(_AtomTable{T})(_atoms(sys, molecule_id = mol.idx))
        @test size(atoms_df(mol), 1) == 0
        @test atoms_df(mol) == atoms_df(sys, molecule_id = mol.idx)
        @test length(atoms(mol)) == 0
        @test atoms(mol) == atoms(sys, molecule_id = mol.idx)
        @test length(collect(eachatom(mol))) == 0
        @test length(collect(eachatom(mol))) == length(collect(eachatom(sys, molecule_id = mol.idx)))
        @test natoms(mol) == 0
        @test natoms(mol) == natoms(sys, molecule_id = mol.idx)

        @test push!(mol, AtomTuple{T}(1, Elements.H)) === mol
        @test length(collect(_atoms(mol))) == 1
        @test size(atoms_df(mol), 1) == 1
        @test size(atoms_df(mol)) == size(atoms_df(sys, molecule_id = mol.idx))
        @test length(atoms(mol)) == 1
        @test atoms(mol) == atoms(sys, molecule_id = mol.idx)
        @test length(collect(eachatom(mol))) == 1
        @test length(collect(eachatom(mol))) == length(collect(eachatom(sys, molecule_id = mol.idx)))
        @test natoms(mol) == 1
        @test natoms(mol) == natoms(sys, molecule_id = mol.idx)

        for atom in eachatom(mol)
            @test parent_molecule(atom) === mol
        end
        @test parent_molecule(Atom(mol, 2, Elements.C)) === mol

        # molecule bonds
        @test size(collect(_bonds(mol)), 1) == 0
        @test Tables.materializer(_BondTable)(_bonds(mol)) ==
            Tables.materializer(_BondTable)(_bonds(sys, molecule_id = mol.idx))
        @test size(bonds_df(mol), 1) == 0
        @test bonds_df(mol) == bonds_df(sys, molecule_id = mol.idx)
        @test length(bonds(mol)) == 0
        @test bonds(mol) == bonds(sys, molecule_id = mol.idx)
        @test length(collect(eachbond(mol))) == 0
        @test length(collect(eachbond(mol))) == length(collect(eachbond(sys, molecule_id = mol.idx)))
        @test nbonds(mol) == 0
        @test nbonds(mol) == nbonds(sys, molecule_id = mol.idx)

        @test push!(mol, BondTuple(
            Atom(mol, 1, Elements.H).idx,
            Atom(mol, 2, Elements.C).idx,
            BondOrder.Single)
        ) === mol
        @test size(collect(_bonds(mol)), 1) == 1
        @test size(collect(_bonds(mol))) == size(collect(_bonds(sys, molecule_id = mol.idx)))
        @test size(bonds_df(mol), 1) == 1
        @test size(bonds_df(mol)) == size(bonds_df(sys, molecule_id = mol.idx))
        @test length(bonds(mol)) == 1
        @test bonds(mol) == bonds(sys, molecule_id = mol.idx)
        @test length(collect(eachbond(mol))) == 1
        @test length(collect(eachbond(mol))) == length(collect(eachbond(sys, molecule_id = mol.idx)))
        @test nbonds(mol) == 1
        @test nbonds(mol) == nbonds(sys, molecule_id = mol.idx)
    end
end
