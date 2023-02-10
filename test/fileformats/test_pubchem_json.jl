@testset "PubChem" begin
        sys = load_pubchem_json("data/aspirin_pug.json")

        @test sys isa System{Float32}
        @test nmolecules(sys) == 1

        mol = molecules(sys)[1]
        @test mol isa Molecule{Float32}
        @test mol.name == "data/aspirin_pug.json_2244"
        @test length(mol.properties) == 3

        @test natoms(mol) == 21
        @test nbonds(mol) == 21


        # used for testing bond orders / file manually manipulated
        sys2 = load_pubchem_json("data/aspirin_pug_bonds.json")

        @test natoms(sys2) == 21
        @test nbonds(sys2) == 21
        @test length(sys2.properties) == 0

        # used for testing bond annotations as properties of bonds
        values = ["BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",       
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC",     
                  "BA_RESONANCE", "BA_BOLD", "BA_FISCHER", "BA_CLOSECONTACT", 
                  "BA_UNKNOWN", "BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC"]
     
        for (i, bond) in enumerate(eachbond(sys2))
                @test bond.properties["PCBondAnnotation_for_conformer"][1] == values[i]
                @test bond.properties["PCBondAnnotation_for_conformer"][2] == values[i]
        end
end
