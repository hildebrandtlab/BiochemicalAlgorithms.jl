@testset "Simple Protein" begin

    mol = Protein()

    @test mol isa Protein{Float32}
    @test mol isa AbstractMolecule{Float32}
    @test name(mol) == ""
    @test atoms(mol) isa SubDataFrame
    @test size(atoms(mol)) == (0, 10)
    @test bonds(mol) isa DataFrame
    @test size(bonds(mol)) == (0, 4)
    @test residues(mol) isa DataFrame
    @test size(residues(mol), 1) == 0
    @test count_atoms(mol) == 0
    @test count_bonds(mol) == 0
    @test length(properties(mol)) == 0
    set_name!(mol, "my_fancy_protein")
    @test name(mol) == "my_fancy_protein"
end

@testset "Filled Protein" begin
    mol = Protein("my_fancy_molecule")

    # add atoms
    r_tmp = Vector3{Float32}(1.0, 2.0, 4.0)
    atoms_df = DataFrame(number = [i for i in 1:6],
                      name = fill("H", 6),
                      element = fill(Elements.H, 6),
                      atomtype = fill("na", 6),
                      r = [i .* Vector3{Float32}(1.0, 2.0, 4.0) for i in 1:6],
                      v = fill(Vector3{Float32}(0.0, 0.0, 0.0), 6), 
                      F = fill(Vector3{Float32}(0.0, 0.0, 0.0), 6),
                      has_velocity = fill(false, 6),
                      has_force = fill(false, 6),
                      properties = Properties()
    )
    r1 = (number = 1, type = AminoAcid('A'), chain = "A")::Residue
    r2 = (number = 2, type = AminoAcid('D'), chain = "B")::Residue
    for (i, atom) in enumerate(eachrow(atoms_df))
        push!(mol, i < 4 ? r1 : r2, copy(atom)::Atom{Float32})
    end

    # add bonds
    for i in 1:4
        push!(mol, (a1 = i, a2 = 5, order = BondOrder.Single, properties = Properties())::Bond)
    end

    @test name(mol) == "my_fancy_molecule"
    @test count_bonds(mol) == 4
    @test count_atoms(mol) == 6
    @test nrow(atoms(mol, r1)) == 3
    @test nrow(atoms(mol, r2)) == 3

    # test properties
    props = properties(mol)
    props["molecule computed"] = false
    props["molecule resolution"] = 2.5
    @test length(props) == 2
    @test length(properties(mol)) == 2

    props["resolution unit"] = "Angstroem"
    @test length(properties(mol)) == 3
    @test !properties(mol)["molecule computed"]
    @test haskey(properties(mol), "resolution unit")
end
