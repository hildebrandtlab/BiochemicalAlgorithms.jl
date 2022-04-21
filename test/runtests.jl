using BALL
using Test

function test_pubchem()
    @testset "PubChem" begin
        mol = load_pubchem_json("./data/aspirin_pug.json")

        @test mol.name == "./data/aspirin_pug.json"

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21
    end
end

@testset "BALL.jl" begin
    test_pubchem()
end
