@testset "Molecule" begin

    # base construction
    mol = Molecule()

    @test mol isa Molecule
    @test mol.name == ""
    @test atoms_df(mol) isa AbstractDataFrame
    @test size(atoms_df(mol)) == (0,14)
    @test length(atoms(mol)) == 0
    @test bonds_df(mol) isa AbstractDataFrame
    @test size(bonds_df(mol)) == (0,5)
    @test length(bonds(mol)) == 0
    @test natoms(mol) == 0
    @test nbonds(mol) == 0
    @test length(mol.properties) == 0
    # set name
    mol.name = "my_fancy_molecule"
    @test mol.name == "my_fancy_molecule"
    
    # add atoms
    for i in 1:6
        atom = ( 
            idx = 0,
            number = 1,
            element = Elements.H,
            name = "my fancy atom",
            atomtype = "heavy",
            r = Vector3{Float32}(i*1.0, i*2.0, i*4.0),
            v = Vector3{Float32}(1.0, 1.0, 1.0),
            F = Vector3{Float32}(0.0, 0.0, 0.0),
            formal_charge = 1,
            charge = 2.0f32,
            radius = 1.02f32,
            has_velocity = true,
            has_force = false,
            properties = Properties()
        )::AtomTuple{Float32}
        push!(mol, atom)
        @test natoms(mol) == i
    end

    # add bonds
    for i in 1:4
        bond = (
            idx = 0,
            a1 = i, 
            a2 = 3, 
            order = BondOrder.Single, 
            properties = Properties()
        )::BondTuple
        push!(mol, bond)
        @test nbonds(mol) == i
    end

    # add properties
    props = mol.properties
    props["molecule computed"] = false 
    props["molecule resolution"] = 2.5
    @test length(props) == 2

    props["resolution unit"] = "Angstroem"
    @test length(props) == 3
    @test !props["molecule computed"]
    @test haskey(props, "resolution unit")
    @test !has_property(mol, "ABCDEFG")

    @test has_property(mol, "resolution unit")
    @test get_property(mol, "resolution unit") == "Angstroem"

    set_property(mol, "test property", "test value")
    @test get_property(mol, "test property") == "test value"
end
