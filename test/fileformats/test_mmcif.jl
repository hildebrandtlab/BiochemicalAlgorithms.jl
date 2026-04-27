@testitem "Read PDBx/mmCIF" begin
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/5pti.cif"), T)
        @test sys isa System{T}
        @test sys.name == "5PTI"  # from data block name
        @test natoms(sys) == 1087
        @test nbonds(sys) == 3  # 3 disulfide bonds from _struct_conn
        @test nmolecules(sys) == 1
        @test nchains(sys) == 2  # ATOM and HETATM split into separate chains, both named "A"
        @test nfragments(sys) == 123
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 58

        # verify disulfide bonds have correct flags
        for b in bonds(sys)
            @test has_flag(b, :TYPE__DISULPHIDE_BOND)
        end

        # verify coordinates of first atom (ARG-1 N)
        a1 = first(atoms(sys))
        @test a1.name == "N"
        @test a1.element == Elements.N
        @test isapprox(a1.r, Vector3{T}(32.231, 15.281, -13.143); atol=T(0.001))

        # secondary structures (2 helices + 3 sheet ranges + coils)
        @test nsecondary_structures(sys) > 0

        # IO loading
        sys2 = open(io -> load_mmcif(io, T), ball_data_path("../test/data/5pti.cif"))
        @test sys2 isa System{T}
        @test sys2.name == "5PTI"  # name comes from data block in both cases
        @test sys.name == sys2.name
        @test natoms(sys2) == 1087
        @test nbonds(sys2) == 3
        @test nmolecules(sys2) == 1
        @test nchains(sys2) == 2
        @test nfragments(sys2) == 123
        @test nnucleotides(sys2) == 0
        @test nresidues(sys2) == 58
    end
end

@testitem "Read PDBx/mmCIF nonexistent" begin
    @test_throws SystemError load_mmcif("nonexistent_file.cif")
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

        # coordinate round-trip
        for (a1, a2) in zip(atoms(sys), atoms(sys2))
            @test isapprox(a1.r, a2.r; atol=T(0.001))
        end

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

@testitem "mmCIF vs PDB equivalence: 5PTI" begin
    # 5PTI is a single-chain hydrolase inhibitor. Expected format differences we
    # don't assert across: sys.name (PDB TITLE vs CIF data block), nbonds (PDB
    # also has CONECT records; CIF has only _struct_conn SS bonds).
    for T in [Float32, Float64]
        sys_pdb = load_pdb(ball_data_path("../test/data/5PTI.pdb"), T)
        sys_cif = load_mmcif(ball_data_path("../test/data/5pti.cif"), T)

        @test natoms(sys_pdb) == natoms(sys_cif)
        @test nchains(sys_pdb) == nchains(sys_cif)
        @test nfragments(sys_pdb) == nfragments(sys_cif)
        @test nfragments.(chains(sys_pdb)) == nfragments.(chains(sys_cif))
        @test nresidues(sys_pdb) == nresidues(sys_cif)
        @test nnucleotides(sys_pdb) == nnucleotides(sys_cif)

        # HETATM flag agreement
        hetero_pdb = count(a -> has_property(a, :is_hetero_atom) &&
                                Bool(get_property(a, :is_hetero_atom)), atoms(sys_pdb))
        hetero_cif = count(a -> has_property(a, :is_hetero_atom) &&
                                Bool(get_property(a, :is_hetero_atom)), atoms(sys_cif))
        @test hetero_pdb == hetero_cif

        # first-atom coordinate agreement (order is deterministic in both parsers)
        a_pdb = first(atoms(sys_pdb))
        a_cif = first(atoms(sys_cif))
        @test a_pdb.name == a_cif.name
        @test a_pdb.element == a_cif.element
        @test isapprox(a_pdb.r, a_cif.r; atol=T(0.001))

        # 3 disulfide bonds should be present in the CIF side
        ss_cif = count(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(sys_cif))
        @test ss_cif == 3
    end
end

@testitem "mmCIF vs PDB equivalence: 2ptc" begin
    # 2ptc is a 4-chain proteinase/inhibitor complex. This entry has no CONECT
    # records in the PDB, so even nbonds agrees between formats (9 disulfides).
    for T in [Float32, Float64]
        sys_pdb = load_pdb(ball_data_path("../test/data/2ptc.pdb"), T)
        sys_cif = load_mmcif(ball_data_path("../test/data/2ptc.cif"), T)

        @test natoms(sys_pdb) == natoms(sys_cif) == 2241
        @test nchains(sys_pdb) == nchains(sys_cif) == 4
        @test nfragments(sys_pdb) == nfragments(sys_cif) == 439
        @test nfragments.(chains(sys_pdb)) == nfragments.(chains(sys_cif)) == [223, 58, 123, 35]
        @test nresidues(sys_pdb) == nresidues(sys_cif)

        # 9 disulfides (6 trypsin + 3 BPTI)
        ss_pdb = count(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(sys_pdb))
        ss_cif = count(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(sys_cif))
        @test ss_pdb == ss_cif == 9

        # HETATM flag agreement
        hetero_pdb = count(a -> has_property(a, :is_hetero_atom) &&
                                Bool(get_property(a, :is_hetero_atom)), atoms(sys_pdb))
        hetero_cif = count(a -> has_property(a, :is_hetero_atom) &&
                                Bool(get_property(a, :is_hetero_atom)), atoms(sys_cif))
        @test hetero_pdb == hetero_cif

        # first-atom coordinate agreement
        a_pdb = first(atoms(sys_pdb))
        a_cif = first(atoms(sys_cif))
        @test a_pdb.name == a_cif.name == "N"
        @test isapprox(a_pdb.r, Vector3{T}(18.871, 65.715, 12.731); atol=T(0.001))
        @test isapprox(a_pdb.r, a_cif.r; atol=T(0.001))
    end
end

@testitem "mmCIF read: 4hhb (HETATMs)" begin
    # 4hhb: hemoglobin tetramer with 4 hemes and ~220 waters. Tests HETATM
    # classification, heme residue recognition, and element parsing for Fe.
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/4hhb.cif"), T)
        @test sys isa System{T}
        @test natoms(sys) == 4779
        @test nchains(sys) == 12
        @test nfragments(sys) == 801

        # hetero-atom count agrees with the PDB parser's result (test_pdb.jl)
        hetero_count = count(
            a -> has_property(a, :is_hetero_atom) && Bool(get_property(a, :is_hetero_atom)),
            atoms(sys))
        @test hetero_count == 395

        # four HEM fragments (one per hemoglobin chain) and four Fe atoms
        hem_frags = filter(f -> f.name == "HEM", collect(fragments(sys)))
        @test length(hem_frags) == 4

        fe_atoms = filter(a -> a.element == Elements.Fe, collect(atoms(sys)))
        @test length(fe_atoms) == 4
        for fe in fe_atoms
            @test fe.name == "FE"
            @test has_property(fe, :is_hetero_atom)
            @test Bool(get_property(fe, :is_hetero_atom))
        end
        # coordinates of the first Fe atom
        @test isapprox(first(fe_atoms).r, Vector3{T}(18.362, 18.488, 23.755); atol=T(0.001))

        # waters (HOH) should all be flagged hetero
        water_frags = filter(f -> f.name == "HOH", collect(fragments(sys)))
        @test length(water_frags) == 221
        for f in water_frags, a in atoms(f)
            @test has_property(a, :is_hetero_atom)
            @test Bool(get_property(a, :is_hetero_atom))
        end
    end
end

@testitem "mmCIF read: 1bna (nucleic acids)" begin
    # 1bna: B-DNA dodecamer (Drew-Dickerson). Tests nucleotide classification
    # and non-protein fragment handling. Note: for this entry the mmCIF parser
    # classifies DA/DT/DG/DC as nucleotides while the PDB parser classifies them
    # as residues — we do not assert residue/nucleotide counts across formats.
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/1bna.cif"), T)
        @test sys isa System{T}
        @test natoms(sys) == 566
        @test nnucleotides(sys) == 24
        @test nresidues(sys) == 0

        # every non-water fragment is a DNA nucleotide
        dna_names = Set(("DA", "DT", "DG", "DC"))
        nuc_frags = filter(f -> f.name in dna_names, collect(fragments(sys)))
        @test length(nuc_frags) == 24

        # the two DNA chains each contain the full 12 nucleotides
        @test nnucleotides.(chains(sys))[1:2] == [12, 12]

        # cross-check against the PDB parser on atom count (both formats agree here)
        sys_pdb = load_pdb(ball_data_path("../test/data/1bna.pdb"), T)
        @test natoms(sys_pdb) == natoms(sys)

        # spot-check: first atom is O5' of the first nucleotide
        a1 = first(atoms(sys))
        @test a1.name == "O5'"
        @test a1.element == Elements.O
        @test isapprox(a1.r, Vector3{T}(18.935, 34.195, 25.617); atol=T(0.001))
    end
end

@testitem "mmCIF real-structure round-trip" begin
    # Round-trip a real multi-fragment structure: 5pti.cif has 58 residues,
    # 3 disulfide bonds, and HETATMs (waters + ligand). The writer currently
    # preserves atoms, chains, fragments, coordinates, disulfide bonds, and
    # HETATM flags — but not CONECT-style bonds or arbitrary metadata.
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/5pti.cif"), T)

        tmpfile = tempname() * ".cif"
        try
            write_mmcif(tmpfile, sys)
            sys2 = load_mmcif(tmpfile, T)

            @test natoms(sys2) == natoms(sys)
            @test nchains(sys2) == nchains(sys)
            @test nfragments(sys2) == nfragments(sys)

            # coordinates survive writer's 3-decimal-place precision
            for (a1, a2) in zip(atoms(sys), atoms(sys2))
                @test isapprox(a1.r, a2.r; atol=T(0.001))
            end

            # disulfide bonds survive
            ss_before = count(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(sys))
            ss_after  = count(b -> has_flag(b, :TYPE__DISULPHIDE_BOND), bonds(sys2))
            @test ss_before == ss_after == 3

            # HETATM flags survive
            hetero_before = count(a -> has_property(a, :is_hetero_atom) &&
                                       Bool(get_property(a, :is_hetero_atom)), atoms(sys))
            hetero_after  = count(a -> has_property(a, :is_hetero_atom) &&
                                       Bool(get_property(a, :is_hetero_atom)), atoms(sys2))
            @test hetero_before == hetero_after
        finally
            rm(tmpfile; force=true)
        end
    end
end
