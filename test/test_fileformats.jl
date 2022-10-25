@testset "PDB" begin
        bpti = load_pdb("data/bpti.pdb")

        @test bpti.name == "bpti.pdb"

        @test count_atoms(bpti) == 454
        @test count_bonds(bpti) == 0

        pdb_5pti = load_pdb("data/5PTI.pdb")

        @test pdb_5pti.name == "5PTI.pdb"

        @test count_atoms(pdb_5pti) == 1087

        @test count_bonds(pdb_5pti) == 0
end

@testset "PubChem" begin
        mol = load_pubchem_json("data/aspirin_pug.json")
        
        @test mol.name == "data/aspirin_pug.json"

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21


        # used for testing bond orders / file manually manipulated
        mol2 = load_pubchem_json("data/aspirin_pug_bonds.json")

        @test count_atoms(mol2) == 21
        @test count_bonds(mol2) == 21
    end

@testset "mol2_import" begin
        mol = load_mol2("data/sustiva.mol2")
        
        @test mol.name == "data/sustiva.mol2"

        @test count_atoms(mol) == 30
        @test count_bonds(mol) == 32
end

@testset "mol2_export" begin

        mol = load_pubchem_json("data/Sustiva_Efavirenz_Conformer3D_CID_64139.json")
        export_mol2(mol)

        exported_mol2_file = open("data/Sustiva_Efavirenz_Conformer3D_CID_64139.mol2")
        @test exported_mol2_file
end