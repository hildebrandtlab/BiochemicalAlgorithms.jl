using BALL
using Test

function test_pubchem()
    @testset "PubChem" begin
        sys = System()
        
        load_pubchem_json!(sys, "./data/aspirin_pug.json")

        @test size(sys.molecules)[1] == 1
        @test sys.molecules[1, :name] == "./data/aspirin_pug.json"

        @test size(sys.atoms)[1] == 21
    end
end

@testset "BALL.jl" begin
    test_pubchem()
end
