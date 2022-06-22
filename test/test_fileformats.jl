@testset "PDB" begin
        bpti = load_pdb("test/data/bpti.pdb")

        @test bpti.name == "bpti.pdb"

        @test count_atoms(bpti) == 454
        @test count_bonds(bpti) == 0

        pdb_5pti = load_pdb("test/data/5PTI.pdb")

        @test pdb_5pti.name == "5PTI.pdb"

        @test count_atoms(pdb_5pti) == 1087

        @test count_bonds(pdb_5pti) == 0
end

@testset "PubChem" begin
        mol = load_pubchem_json("test/data/aspirin_pug.json")
        
        @test mol.name == "test/data/aspirin_pug.json"

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21
    end