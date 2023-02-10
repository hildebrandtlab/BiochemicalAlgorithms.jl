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