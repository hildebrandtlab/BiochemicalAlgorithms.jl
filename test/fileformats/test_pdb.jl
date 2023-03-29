@testitem "PDB" begin
    bpti = load_pdb(ball_data_path("../test/data/bpti.pdb"))

    @test bpti.name == "bpti.pdb"

    @test natoms(bpti) == 454
    @test nbonds(bpti) == 0

    pdb_5pti = load_pdb(ball_data_path("../test/data/5PTI.pdb"))

    @test pdb_5pti.name == "5PTI.pdb"

    @test natoms(pdb_5pti) == 1087

    @test nbonds(pdb_5pti) == 0
end
