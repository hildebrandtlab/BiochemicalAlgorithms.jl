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
        @test mol.name == "data/aspirin_pug.json_2244"
        @test length(mol.properties) == 3

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21


        # used for testing bond orders / file manually manipulated
        mol2 = load_pubchem_json("data/aspirin_pug_bonds.json")[1]

        @test count_atoms(mol2) == 21
        @test count_bonds(mol2) == 21
        @test length(mol2.properties) == 0

        # used for testing bond annotations as properties of bonds
        values = ["BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",       
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC",     
                  "BA_RESONANCE", "BA_BOLD", "BA_FISCHER", "BA_CLOSECONTACT", 
                  "BA_UNKNOWN", "BA_CROSSED", "BA_DASHED","BA_WAVY", "BA_DOTTED",
                  "BA_WEDGE_UP", "BA_WEDGE_DOWN", "BA_ARROW", "BA_AROMATIC"]
     
        for i in eachindex(mol2.bonds.properties) 
                @test mol2.bonds.properties[i]["PCBondAnnotation_for_conformer"][1] == values[i]
                @test mol2.bonds.properties[i]["PCBondAnnotation_for_conformer"][2] == values[i]
        end
end


@testset "mol2_import" begin
        mol = load_mol2("data/sustiva_openbabel.mol2")
        
        @test mol.name == "data/sustiva_openbabel.mol2"

        @test count_atoms(mol) == 30
        @test count_bonds(mol) == 32
end

@testset "mol2_export" begin
        export_mol2(load_pubchem_json("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.json"), "data/")

        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139_balljl_export.mol2")[1] == "@<TRIPOS>MOLECULE"
        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139_balljl_export.mol2")[8] == "@<TRIPOS>ATOM"
        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139_balljl_export.mol2")[39] == "@<TRIPOS>BOND"

        rm("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139_balljl_export.mol2")

end