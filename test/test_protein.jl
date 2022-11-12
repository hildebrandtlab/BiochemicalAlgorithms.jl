@testset "Simple Protein" begin

    protein = Protein("my_fancy_protein")

    @test protein isa Protein
    @test protein isa AbstractMolecule
    @test !isa(protein, Molecule)
    @test protein.name == "my_fancy_protein"
    @test protein.atoms isa DataFrame
    @test size(protein.atoms) == (0, 14)
    @test protein.bonds isa DataFrame
    @test size(protein.bonds) == (0, 3)
    @test protein.chains isa Vector{ProteinChain}
    @test size(protein.chains, 1) == 0
    @test count_atoms(protein) == 0
    @test count_bonds(protein) == 0
end

@testset "Filled Protein" begin

    protein = Protein()
    # now create PDBChains
    chain_A = ProteinChain("chain A", 
                            DataFrame(number = [1, 2, 3],
                                      type = fill(AminoAcid(0), 3),
                                      chain_id = ["A", "A", "A"] ))

    chain_B = ProteinChain("chain B", 
                            DataFrame(number = [1, 2, 3],
                                      type = fill(AminoAcid(3), 3),
                                      chain_id = ["B", "B", "B"] ))
                                      
    push!(protein.chains, chain_A)
    push!(protein.chains, chain_B)
    @test size(protein.chains, 1) == 2
    @test protein.chains[1] isa ProteinChain
end
