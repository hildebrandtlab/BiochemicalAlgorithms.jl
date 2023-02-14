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
        mols = load_pubchem_json("data/aspirin_pug.json")

        @test mols isa Vector{Molecule}
        @test length(mols) == 1

        mol = mols[1]
        @test mol isa Molecule
        @test name(mol) == "data/aspirin_pug.json_2244"
        @test length(properties(mol)) == 3

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21


        # used for testing bond orders / file manually manipulated
        mol2 = load_pubchem_json("data/aspirin_pug_bonds.json")[1]

        @test count_atoms(mol2) == 21
        @test count_bonds(mol2) == 21
        @test length(properties(mol2)) == 0

        # used for testing bond annotations as properties of bonds
        values = ["BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",       
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC",     
                  "BA_RESONANCE", "BA_BOLD", "BA_FISCHER", "BA_CLOSECONTACT", 
                  "BA_UNKNOWN", "BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC"]
     
        for i in eachindex(bonds(mol2).properties) 
                @test bonds(mol2).properties[i]["PCBondAnnotation_for_conformer"][1] == values[i]
                @test bonds(mol2).properties[i]["PCBondAnnotation_for_conformer"][2] == values[i]
        end
end
