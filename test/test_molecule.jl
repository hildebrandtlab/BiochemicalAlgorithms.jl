@testset "Molecule" begin

    # base construction
    mol = Molecule()

    @test mol isa Molecule
    @test mol.name == ""
    @test mol.atoms isa DataFrame
    @test size(mol.atoms) == (0,11)
    @test mol.bonds isa DataFrame
    @test size(mol.bonds) == (0,3)
    @test count_atoms(mol) == 0
    @test count_bonds(mol) == 0
    # set name
    mol.name = "my_fancy_molecule"
    @test mol.name == "my_fancy_molecule"
    
    # add atoms
    for i in 1:6
        atom = ( number = 1,
        name = "my fancy atom", 
        element = Elements.H, 
        atomtype = "heavy",
        r = Vector3{Float64}(i*1.0, i*2.0, i*4.0),
        v = Vector3{Float64}(1.0, 1.0, 1.0),
        F = Vector3{Float64}(0.0, 0.0, 0.0),
        has_velocity = true,
        has_force = false,
        frame_id = 1,
        properties = Dict{String, Any}()
        )
        push!(mol, atom)
        @test count_atoms(mol) == i
    end

    # add bonds
    for i in 1:4
        bond = (a1 = i, 
                a2 = 3, 
                order = BondOrder.Single)
        push!(mol, bond)
        @test count_bonds(mol) == i
    end
end
