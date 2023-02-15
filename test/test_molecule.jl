@testset "Molecule" begin

    # base construction
    mol = Molecule()

    @test mol isa Molecule
    @test name(mol) == ""
    @test atoms(mol) isa SubDataFrame
    @test size(atoms(mol)) == (0,10)
    @test bonds(mol) isa DataFrame
    @test size(bonds(mol)) == (0,4)
    @test count_atoms(mol) == 0
    @test count_bonds(mol) == 0
    @test length(properties(mol)) == 0
    # set name
    set_name!(mol, "my_fancy_molecule")
    @test name(mol) == "my_fancy_molecule"
    
    # add atoms
    for i in 1:6
        atom = ( number = 1,
            name = "my fancy atom", 
            element = Elements.H, 
            atomtype = "heavy",
            r = Vector3{Float32}(i*1.0, i*2.0, i*4.0),
            v = Vector3{Float32}(1.0, 1.0, 1.0),
            F = Vector3{Float32}(0.0, 0.0, 0.0),
            has_velocity = true,
            has_force = false,
            properties = Properties()
        )::Atom{Float32}
        push!(mol, atom)
        @test count_atoms(mol) == i
    end

    # add bonds
    for i in 1:4
        bond = (a1 = i, 
            a2 = 3, 
            order = BondOrder.Single, 
            properties = Properties()
        )::Bond
        push!(mol, bond)
        @test count_bonds(mol) == i
    end

    # add properties
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
