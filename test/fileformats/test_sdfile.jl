@testset "SDFile" begin
    @testset "Reading" begin
        mols = molecules(load_sdfile("data/sdfile_test_1.sdf"))

        @test length(mols) == 11

        @test sum(map(natoms, mols)) == 518
        @test sum(map(nbonds, mols)) == 528

        mol = mols[1]
        @test natoms(mol) == 39
        @test nbonds(mol) == 42

        @test has_property(mol, "NAME")
        @test get_property(mol, "NAME") == "Abacavir_sulfate"

        @test has_property(mol, "b_1rotN")
        @test get_property(mol, "b_1rotN") == "6"

        @test has_property(mol, "Weight")
        @test get_property(mol, "Weight") == "286.339"

        @test has_property(mol, "TPSA")
        @test get_property(mol, "TPSA") == "101.88"

        @test has_property(mol, "a_acc")
        @test get_property(mol, "a_acc") == "4"

        @test has_property(mol, "a_don")
        @test get_property(mol, "a_don") == "3"

        @test has_property(mol, "logP(o/w)")
        @test get_property(mol, "logP(o/w)") == "0.40906"

        @test has_property(mol, "SlogP")
        @test get_property(mol, "SlogP") == "1.1878"
    
    end
end