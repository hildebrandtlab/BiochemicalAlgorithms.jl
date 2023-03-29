@testitem "PDB" begin
    using BioStructures:
    PDB,
    countatoms,
    countchains,
    countresidues,
    read

    sys = load_pdb(ball_data_path("../test/data/bpti.pdb"))
    @test sys.name == "bpti.pdb"
    @test natoms(sys) == 454
    @test nbonds(sys) == 0
    @test nmolecules(sys) == 1
    @test nchains(sys) == 1
    @test nfragments(sys) == 58
    @test nresidues(sys) == 0
    @test nnucleotides(sys) == 0

    sys = load_pdb(ball_data_path("../test/data/5PTI.pdb"))
    @test sys.name == "5PTI.pdb"
    @test natoms(sys) == 1087
    @test nbonds(sys) == 0
    @test nmolecules(sys) == 1
    @test nchains(sys) == 1
    @test nfragments(sys) == 123
    @test nresidues(sys) == 0
    @test nnucleotides(sys) == 0

    sys = load_pdb("data/2ptc.pdb")
    orig = read("data/2ptc.pdb", PDB)
    orig_model = orig.models[1]
    @test sys.name == orig.name
    @test natoms(sys) == countatoms(orig)
    @test nbonds(sys) == 0
    @test nmolecules(sys) == 1
    @test nchains(sys) == countchains(orig)
    @test nfragments(sys) == countresidues(orig)
    for chain in eachchain(sys)
        orig_chain = get(orig_model.chains, chain.name, nothing)
        @test !isnothing(orig_chain)
        @test natoms(chain) == countatoms(orig_chain)
        @test nfragments(chain) == countresidues(orig_chain)
        for frag in eachfragment(chain)
            # reconstruct BioStructures naming for residues
            orig_name = strip("$(get_property(frag, :is_hetero_fragment, false) ? "H_" : "")\
                $(frag.number)$(get_property(frag, :insertion_code, ""))")
            orig_frag = get(orig_chain.residues, orig_name, nothing)
            @test !isnothing(orig_frag)
            @test natoms(frag) == countatoms(orig_frag)
        end
    end
end
