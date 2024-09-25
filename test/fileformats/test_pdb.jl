@testitem "PDB reading" begin
	s = load_pdb(ball_data_path("../test/data/pdbfile_test.pdb"))

    @test natoms(s) == 892  
    @test nfragments(s) == 58
    @test nresidues(s) == 58
    @test nnucleotides(s) == 0
	@test nchains(s) == 1
    @test nsecondary_structures(s) == 7

    @test [s.type for s in secondary_structures(s)] == [
        SecondaryStructureElement.Coil, SecondaryStructureElement.Strand,
        SecondaryStructureElement.Coil, SecondaryStructureElement.Strand,
        SecondaryStructureElement.Coil, SecondaryStructureElement.Helix,
        SecondaryStructureElement.Coil]


    # test handling of invalid PDB lines
    @test_throws Exception load_pdb(ball_data_path("../test/data/bpti.pdb"); strict_line_checking=true)

    # test auto-correct of invalid PDB lines
    sys = load_pdb(ball_data_path("../test/data/bpti.pdb"); strict_line_checking=false)
    @test sys.name == "PTI (from 2PTC.BRK)"
    @test natoms(sys) == 454
    @test nbonds(sys) == 3
    @test nmolecules(sys) == 1
    @test nchains(sys) == 1
    @test nfragments(sys) == 58
    @test nnucleotides(sys) == 0
    @test nresidues(sys) == 58

    sys = load_pdb(ball_data_path("../test/data/5PTI.pdb"))
    @test sys.name == "HYDROLASE INHIBITOR"
    @test natoms(sys) == 1087
    @test nbonds(sys) == 7
    @test nmolecules(sys) == 1
    @test nchains(sys) == 2
    @test nfragments(sys) == 123
    @test nnucleotides(sys) == 0
    @test nresidues(sys) == 58

    sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
    @test natoms(sys) == 2241
    @test nbonds(sys) == 9
    @test nmolecules(sys) == 1
    @test nchains(sys) == 4
    @test nfragments(sys) == 439
    @test nfragments.(chains(sys)) == [223, 58, 123, 35]
    @test nnucleotides.(chains(sys)) == [0, 0, 0, 0]
    @test nresidues.(chains(sys)) == [223, 58, 0, 0] 



    
end
