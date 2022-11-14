@testset "Simple PDBMolecule" begin

    mol = PDBMolecule("my_fancy_molecule",)

    @test mol isa PDBMolecule
    @test mol isa AbstractMolecule
    @test !isa(mol, Molecule)
    @test mol.name == "my_fancy_molecule"
    @test mol.atoms isa DataFrame
    @test size(mol.atoms) == (0, 14)
    @test mol.bonds isa DataFrame
    @test size(mol.bonds) == (0, 4)
    @test mol.chains isa Vector{PDBChain}
    @test size(mol.chains, 1) == 0
    @test count_atoms(mol) == 0
    @test count_bonds(mol) == 0
end  

@testset "Filled PDBMolecule" begin

    bonds = DataFrame(a1 = [1, 2, 3, 4], a2 = fill(3, 4), order = fill(BondOrder.Single, 4))

    r_tmp = Vector3{Float64}(1.0, 2.0, 4.0)
 
    atoms = DataFrame(number = [i for i in 1:6],
                      name = fill("H", 6),
                      element = fill(Elements.H, 6),
                      atomtype = fill("na", 6),
                      r = [i * Vector3{Float64}(1.0, 2.0, 4.0) for i in 1:6],
                      v = fill(Vector3{Float64}(0.0, 0.0, 0.0), 6), 
                      F = fill(Vector3{Float64}(0.0, 0.0, 0.0), 6),
                      has_velocity = fill(false, 6),
                      has_force = fill(false, 6),
                      frame_id  = fill(1, 6))
    
    mol = PDBMolecule("my_fancy_molecule", atoms, bonds)
    
    @test mol.name == "my_fancy_molecule"
    @test count_bonds(mol) == 4
    @test count_atoms(mol) == 6


    # now create PDBChains
    chain_A = PDBChain("chain A", DataFrame(number = [1,2,3],
                                            name = ["f1", "f2", "f3"],
                                            chain_id = ["A", "A", "A"] ))
    chain_B = PDBChain("chain B", DataFrame(number = [1,2,3],
                                            name = ["f1", "f2", "f3"],
                                            chain_id = ["B", "B", "B"] ))
    push!(mol.chains, chain_A)
    push!(mol.chains, chain_B)
    @test size(mol.chains, 1) == 2
    @test mol.chains[1] isa PDBChain

end

