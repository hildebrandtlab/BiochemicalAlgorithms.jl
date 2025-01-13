@testitem "Read PDB" begin
    using BioStructures:
        PDBFormat,
        countatoms,
        countchains,
        countresidues,
        read

    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        @test sys isa System{T}
        @test sys.name == "bpti.pdb"
        @test natoms(sys) == 454
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 58
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58

        sys = load_pdb(ball_data_path("../test/data/5PTI.pdb"), T)
        @test sys isa System{T}
        @test sys.name == "5PTI.pdb"
        @test natoms(sys) == 1087
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 123
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58

        sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"), T)
        orig = read(ball_data_path("../test/data/2ptc.pdb"), PDBFormat)
        orig_model = orig.models[1]
        @test sys isa System{T}
        @test sys.name == orig.name
        @test natoms(sys) == countatoms(orig)
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 1
        @test nchains(sys) == countchains(orig)
        @test nfragments(sys) == countresidues(orig)
        @test nfragments.(chains(sys)) == [346, 93]
        @test nnucleotides.(chains(sys)) == [0, 0]
        @test nresidues.(chains(sys)) == [223, 58]
        for chain in chains(sys)
            orig_chain = get(orig_model.chains, chain.name, nothing)
            @test !isnothing(orig_chain)
            @test natoms(chain) == countatoms(orig_chain)
            @test nfragments(chain) == countresidues(orig_chain)
            for frag in fragments(chain)
                # reconstruct BioStructures naming for residues
                orig_name = strip("$(get_property(frag, :is_hetero_fragment, false) ? "H_" : "")\
                    $(frag.number)$(get_property(frag, :insertion_code, ""))")
                orig_frag = get(orig_chain.residues, orig_name, nothing)
                @test !isnothing(orig_frag)
                @test natoms(frag) == countatoms(orig_frag)
            end
        end
    end
end

@testitem "Write PDB" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol; name = "A")
        frag = Fragment(chain, 1)
        Atom(frag, 1, Elements.C; name = "C")
        Atom(frag, 2, Elements.O; name = "O", r = ones(Vector3{T}))

        (fname, fh) = mktemp(; cleanup = true)

        write_pdb(fname, sys)
        sys2 = load_pdb(fname, T)
        @test sys2 isa System{T}
        @test natoms(sys2) == 2
        @test atoms(sys2).name == ["C", "O"]
        @test nfragments(sys2) == 1
        @test nchains(sys2) == 1
        @test first(chains(sys2)).name == "A"
        @test nmolecules(sys2) == 1
    end
end

@testitem "Read PDBx/mmCIF" begin
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/5pti.cif"), T)
        @test sys isa System{T}
        @test sys.name == "5pti.cif"
        @test natoms(sys) == 1087
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 123
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58
    end
end

@testitem "Write PDBx/mmCIF" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol; name = "A")
        frag = Fragment(chain, 1; name = "FRG1")
        Atom(frag, 1, Elements.C; name = "C")
        Atom(frag, 2, Elements.O; name = "O", r = ones(Vector3{T}))

        (fname, fh) = mktemp(; cleanup = true)

        write_mmcif(fname, sys)
        sys2 = load_mmcif(fname, T)
        @test sys2 isa System{T}
        @test natoms(sys2) == 2
        @test atoms(sys2).name == ["C", "O"]
        @test nfragments(sys2) == 1
        @test nchains(sys2) == 1
        @test first(chains(sys2)).name == "A"
        @test nmolecules(sys2) == 1
    end
end
