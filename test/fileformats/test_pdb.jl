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
        @test sys.name == "PTI (from 2PTC.BRK)"
        @test natoms(sys) == 454
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 58
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58

        sys = load_pdb(ball_data_path("../test/data/5PTI.pdb"), T)
        @test sys isa System{T}
        @test sys.name == "HYDROLASE INHIBITOR"
        @test natoms(sys) == 1087
        @test nbonds(sys) == 7
        @test nmolecules(sys) == 1
        @test nchains(sys) == 2
        @test nfragments(sys) == 123
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58

        sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"), T)
        @test sys isa System{T}
        @test sys.name == "COMPLEX (PROTEINASE/INHIBITOR)"
        @test natoms(sys) == 2241
        @test nbonds(sys) == 9
        @test nmolecules(sys) == 1
        @test nchains(sys) == 4
        @test nfragments(sys) == 439
        @test nfragments.(chains(sys)) == [223, 58, 123, 35]
        @test nnucleotides.(chains(sys)) == [0, 0, 0, 0]
        @test nresidues.(chains(sys)) == [223, 58, 0, 0]
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
