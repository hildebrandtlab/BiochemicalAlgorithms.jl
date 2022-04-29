using BiochemicalAlgorithms
using Test

function test_pubchem()
    @testset "PubChem" begin
        mol = load_pubchem_json("./data/aspirin_pug.json")

        @test mol.name == "./data/aspirin_pug.json"

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21
    end
end

function test_pdb()
    @testset "PDB" begin
        bpti = load_pdb("./data/bpti.pdb")

        @test bpti.name == "bpti.pdb"
        
        @test count_atoms(bpti) == 454
        @test count_bonds(bpti) == 0

        pdb_5pti = load_pdb("./data/5PTI.pdb")

        @test pdb_5pti.name == "5PTI.pdb"

        @test count_atoms(pdb_5pti) == 1087 # currently fails b/c auf alternate locations
        @test count_bonds(pdb_5pti) == 0
    end
end

@testset "BiochemicalAlgorithms.jl" begin
    test_pdb()
    test_pubchem()
end
