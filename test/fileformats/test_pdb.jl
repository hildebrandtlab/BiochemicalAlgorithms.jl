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

        sys2 = open(io -> load_pdb(io, T), ball_data_path("../test/data/bpti.pdb"))
        @test sys == sys2

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

        open(io -> write_pdb(io, sys), fname, "w")
        sys2 = load_pdb(fname, T)
        @test sys2 isa System{T}
        @test natoms(sys2) == 2
        @test atoms(sys2).name == ["C", "O"]
        @test nfragments(sys2) == 1
        @test nchains(sys2) == 1
        @test first(chains(sys2)).name == "A"
        @test nmolecules(sys2) == 1
    end


    # test for single chain
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol; name = "A")
        chain2 = Chain(mol; name = "B")
        frag = Fragment(chain, 1)
        frag2 = Fragment(chain2, 2)
        Atom(frag, 1, Elements.C; name = "C")
        Atom(frag, 2, Elements.O; name = "O", r = ones(Vector3{T}))

        Atom(frag2, 1, Elements.C; name = "C")
        Atom(frag2, 2, Elements.O; name = "H", r = ones(Vector3{T}))
        Atom(frag2, 3, Elements.O; name = "H", r = ones(Vector3{T}))
        (fname, fh) = mktemp(; cleanup = true)

        write_pdb(fname, chain2)
        sys2 = load_pdb(fname, T)
        @test sys2 isa System{T}
        @test natoms(sys2) == 3
        @test atoms(sys2).name == ["C", "H", "H"]
        @test nfragments(sys2) == 1
        @test nchains(sys2) == 1
        @test first(chains(sys2)).name == "B"
        @test nmolecules(sys2) == 1
    end

    # test for single fragment (only coordinates)
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol; name = "A")
        chain2 = Chain(mol; name = "B")
        frag = Fragment(chain, 1)
        frag2 = Fragment(chain2, 2)
        Atom(frag, 1, Elements.C; name = "C")
        Atom(frag, 2, Elements.O; name = "O", r = ones(Vector3{T}))

        Atom(frag2, 1, Elements.C; name = "C")
        Atom(frag2, 2, Elements.O; name = "H", r = ones(Vector3{T}))
        Atom(frag2, 3, Elements.O; name = "H", r = ones(Vector3{T}))
        (fname, fh) = mktemp(; cleanup = true)

        write_pdb(fname, chain2)
        sys2 = load_pdb(fname, T)
        @test sys2 isa System{T}
        @test natoms(sys2) == 3
        @test atoms(sys2).name == ["C", "H", "H"]
        @test nfragments(sys2) == 1
        @test nchains(sys2) == 1
        @test first(chains(sys2)).name == "B"
        @test nmolecules(sys2) == 1
    end
end

@testitem "Read PDB/4hhb" begin
    for T in [Float32, Float64]
        # 4hhb: multi-chain hemoglobin (4 protein chains + waters/heteroatoms)
        sys = load_pdb(ball_data_path("../test/data/4hhb.pdb"), T)
        @test sys isa System{T}
        @test natoms(sys) == 4779
        @test nchains(sys) == 12
        @test nfragments(sys) == 801

        # HETATM records flagged
        hetero_count = count(
            a -> has_property(a, :is_hetero_atom) && Bool(get_property(a, :is_hetero_atom)),
            atoms(sys))
        @test hetero_count == 395
    end
end

@testitem "PDB round-trip coordinates" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        tmpfile = tempname() * ".pdb"

        try
            write_pdb(tmpfile, sys)
            sys2 = load_pdb(tmpfile, T)

            @test natoms(sys2) == natoms(sys)
            @test nfragments(sys2) == nfragments(sys)

            # verify coordinate round-trip (PDB has 3 decimal places)
            for (a1, a2) in zip(atoms(sys), atoms(sys2))
                @test isapprox(a1.r, a2.r; atol=T(0.001))
            end
        finally
            rm(tmpfile; force=true)
        end
    end
end

@testitem "Read PDB/nonexistent" begin
    @test_throws SystemError load_pdb("nonexistent_file.pdb")
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

        sys = open(io -> load_mmcif(io, T), ball_data_path("../test/data/5pti.cif"))
        @test sys isa System{T}
        @test_broken sys.name == "5pti.cif" # BioStructures does not set a name here...
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

        open(io -> write_mmcif(io, sys), fname, "w")
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
