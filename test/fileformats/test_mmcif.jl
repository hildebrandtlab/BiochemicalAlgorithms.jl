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

        # verify disulfide bonds have correct flags (also flagged covalent)
        for b in bonds(sys)
            @test has_flag(b, :TYPE__DISULPHIDE_BOND)
            @test has_flag(b, :TYPE__COVALENT)
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

        # CIF additionally parses 6 metal-coordination bonds (Ca²⁺ in trypsin)
        # from _struct_conn (the PDB reader doesn't); they're flagged covalent.
        @test count(b -> has_flag(b, :TYPE__COVALENT), bonds(sys_cif)) == 15

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

        # 4 metal-coordination bonds (Fe — NE2 of His) parsed from _struct_conn
        @test nbonds(sys) == 4
        for b in bonds(sys)
            @test has_flag(b, :TYPE__COVALENT)
            a1, a2 = get_partners(b)
            elements = (a1.element, a2.element)
            @test Elements.Fe in elements
            @test Elements.N in elements
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
    # 3 disulfide bonds, secondary structure (2 helices + 3 strands + coils),
    # and HETATMs (waters + ligand).
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

            # secondary structures survive in count and per-type composition,
            # and the strand names retain the original sheet:range_id format.
            n_helix(s) = count(ss -> ss.type == SecondaryStructureElement.Helix, secondary_structures(s))
            n_strand(s) = count(ss -> ss.type == SecondaryStructureElement.Strand, secondary_structures(s))
            @test n_helix(sys2)  == n_helix(sys)  == 2
            @test n_strand(sys2) == n_strand(sys) == 3

            strand_names(s) = [ss.name for ss in secondary_structures(s) if ss.type == SecondaryStructureElement.Strand]
            @test strand_names(sys2) == strand_names(sys)
        finally
            rm(tmpfile; force=true)
        end
    end
end

@testitem "mmCIF read: 1bna hydrogen bonds" begin
    # 1bna: B-DNA dodecamer. The mmCIF file's `_struct_conn` loop encodes the
    # base-pair hydrogen bonds; verify they're parsed with TYPE__HYDROGEN.
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/1bna.cif"), T)
        @test nbonds(sys) == 32
        for b in bonds(sys)
            @test has_flag(b, :TYPE__HYDROGEN)
            @test !has_flag(b, :TYPE__COVALENT)
        end
    end
end

@testitem "mmCIF read: NMR multi-model" begin
    # 1d3z: NMR ensemble of ubiquitin (10 models, 1231 atoms each). Verify
    # that all models are loaded as separate frames and that fragments are
    # shared across models — i.e. nfragments matches one model's residue
    # count, not 10×.
    for T in [Float32, Float64]
        sys = load_mmcif(ball_data_path("../test/data/1d3z.cif"), T)
        @test sys isa System{T}
        @test get_property(sys, :experimental_method, "") == "SOLUTION NMR"

        frames = sort(unique(sys._atoms.frame_id))
        @test length(frames) == 10
        @test frames == collect(1:10)

        # Each frame contributes the same atom count.
        per_frame = [count(==(f), sys._atoms.frame_id) for f in frames]
        @test all(==(per_frame[1]), per_frame)
        @test sum(per_frame) == length(sys._atoms.idx)

        # Fragments are shared, not replicated per model.
        @test nfragments(sys) == 76
        @test nchains(sys) == 1

        # natoms() returns just the current frame.
        @test natoms(sys) == per_frame[1] == 1231
    end
end

@testitem "mmCIF metadata extraction" begin
    sys = load_mmcif(ball_data_path("../test/data/5pti.cif"))
    @test get_property(sys, :entry_id, "") == "5PTI"
    @test occursin("BOVINE PANCREATIC TRYPSIN INHIBITOR", get_property(sys, :title, ""))
    @test get_property(sys, :keywords, "") == "Hydrolase Inhibitor"
    @test get_property(sys, :deposition_date, "") == "1984-10-05"
    @test get_property(sys, :experimental_method, "") == "X-RAY DIFFRACTION"
    @test get_property(sys, :resolution, 0.0) ≈ 1.0
    @test get_property(sys, :space_group, "") == "P 21 21 21"
    @test get_property(sys, :cell_a, 0.0) ≈ 74.1
    @test get_property(sys, :cell_b, 0.0) ≈ 23.4
    @test get_property(sys, :cell_c, 0.0) ≈ 28.9
    @test get_property(sys, :cell_alpha, 0.0) ≈ 90.0

    # 4hhb has a non-orthorhombic angle, useful as a contrasting check.
    sys = load_mmcif(ball_data_path("../test/data/4hhb.cif"))
    @test get_property(sys, :entry_id, "") == "4HHB"
    @test get_property(sys, :keywords, "") == "OXYGEN TRANSPORT"
    @test get_property(sys, :resolution, 0.0) ≈ 1.74
    @test get_property(sys, :cell_beta, 0.0) ≈ 99.34
end

@testitem "mmCIF parser: CIF 2.0 syntax" begin
    # Smoke-test the CIF 2.0 magic line and triple-quoted string parsing.
    # Triple-quoted text blocks are recognised when the value starts on the
    # line after the tag (i.e. the first character of the value-line is `\"\"\"`).
    cif = """#\\#CIF_2.0
data_test
_entry.id TEST
_struct.title
\"\"\"a multi-line
title with 'apostrophes' inside\"\"\"
"""
    sys = load_mmcif(IOBuffer(cif))
    @test sys.name == "test"
    @test get_property(sys, :entry_id, "") == "TEST"
    title = get_property(sys, :title, "")
    @test occursin("apostrophes", title)
    @test occursin("multi-line", title)
end

@testitem "mmCIF parser: multi-block CIF picks first" begin
    # Multi-block files are valid CIF; the loader picks the first data block
    # (and warns? actually no — silently picks first per spec).
    cif = """data_first
_entry.id FIRST
data_second
_entry.id SECOND
"""
    sys = load_mmcif(IOBuffer(cif))
    @test get_property(sys, :entry_id, "") == "FIRST"
end

@testitem "mmCIF read: insertion codes" begin
    # Synthetic CIF with two residues at sequence number 100, distinguished by
    # insertion code A and B (a common antibody/immunoglobulin pattern).
    cif = """data_INS
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 100 A 0.000 0.000 0.000 1
ATOM 2 C CA . ALA A 100 A 1.000 0.000 0.000 1
ATOM 3 N N . GLY A 100 B 2.000 0.000 0.000 1
ATOM 4 C CA . GLY A 100 B 3.000 0.000 0.000 1
"""
    sys = load_mmcif(IOBuffer(cif))
    @test natoms(sys) == 4
    @test nfragments(sys) == 2
    frags = collect(fragments(sys))
    @test [f.name for f in frags] == ["ALA", "GLY"]
    @test [get_property(f, :insertion_code, "") for f in frags] == ["A", "B"]
    # Same sequence number, different insertion codes — must be treated as
    # distinct residues, not collapsed.
    @test all(f.number == 100 for f in frags)
end

@testitem "mmCIF read: covale bonds" begin
    # Synthetic CIF with a covale (covalent linker) connection.
    cif = """data_COV
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N  . LYS A 1 ? 0.000 0.000 0.000 1
ATOM 2 C CA . LYS A 1 ? 1.000 0.000 0.000 1
ATOM 3 N NZ . LYS A 1 ? 2.000 0.000 0.000 1
ATOM 4 C C  . PYR A 2 ? 3.000 0.000 0.000 1
ATOM 5 O O  . PYR A 2 ? 4.000 0.000 0.000 1
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
covale1 covale A LYS 1 NZ A PYR 2 C doub
"""
    sys = load_mmcif(IOBuffer(cif))
    @test nbonds(sys) == 1
    b = first(bonds(sys))
    @test has_flag(b, :TYPE__COVALENT)
    @test b.order == BondOrder.Double
    @test get_property(b, :CIF_CONN_TYPE, "") == "covale"
end

@testitem "mmCIF malformed: missing data block" begin
    @test_throws ErrorException load_mmcif(IOBuffer("# just a comment\n"))
end

@testitem "mmCIF parser: blank/comment lines preserved inside text blocks" begin
    # CIF semicolon text blocks may contain blank lines and lines starting
    # with `#`; the parser must not silently drop them.
    cif = """data_TXT
_struct.title
;line one

# this isn't a comment, it's part of the value
line three
;
"""
    sys = load_mmcif(IOBuffer(cif))
    title = get_property(sys, :title, "")
    @test occursin("line one", title)
    @test occursin("line three", title)
    @test occursin("# this isn't a comment", title)
    # blank line between "line one" and the `#` line should also survive
    @test occursin("\n\n", title) || occursin("line one\n\n", title)
end

@testitem "mmCIF parser: triple-quoted on a single line" begin
    # CIF 2.0: a triple-quoted string value may both open and close on the
    # same line. Tested with the tag-then-value-on-next-line form, which is
    # what triggers the parser's `ExpectValue` triple-quote branch.
    cif = """#\\#CIF_2.0
data_TQ
_entry.id TQ
_struct.title
\"\"\"single-line title\"\"\"
"""
    sys = load_mmcif(IOBuffer(cif))
    @test get_property(sys, :entry_id, "") == "TQ"
    @test get_property(sys, :title, "") == "single-line title"

    # Empty triple-quoted body is also valid.
    cif_empty = """#\\#CIF_2.0
data_TQ2
_entry.id TQ2
_struct.title
\"\"\"\"\"\"
"""
    sys = load_mmcif(IOBuffer(cif_empty))
    @test get_property(sys, :entry_id, "") == "TQ2"
    @test get_property(sys, :title, "") == ""
end

@testitem "mmCIF writer: Elements.Unknown maps to ?" begin
    sys = System{Float32}()
    mol = Molecule(sys)
    chain = Chain(mol; name = "A")
    frag = Fragment(chain, 1; name = "UNK")
    Atom(frag, 1, Elements.Unknown; name = "X")

    tmp = tempname() * ".cif"
    try
        write_mmcif(tmp, sys)
        text = read(tmp, String)
        # The ATOM row's third whitespace-separated field is type_symbol.
        atom_lines = filter(l -> startswith(l, "ATOM ") || startswith(l, "HETATM "), split(text, '\n'))
        @test length(atom_lines) == 1
        type_sym = split(atom_lines[1])[3]
        @test type_sym == "?"

        # Round-trip should still load (the reader's parse_element_string
        # accepts `?` and treats it as Unknown).
        sys2 = load_mmcif(tmp)
        @test natoms(sys2) == 1
        @test first(atoms(sys2)).element == Elements.Unknown
    finally
        rm(tmp; force=true)
    end
end

@testitem "mmCIF writer: multiline string in property is collapsed" begin
    # If a property contains a newline, the writer must keep the row on a
    # single line (semicolon text blocks can't be embedded in a row that's
    # written via space-joined fields).
    sys = System{Float32}()
    mol = Molecule(sys)
    chain = Chain(mol; name = "A")
    frag = Fragment(chain, 1; name = "ALA")
    Atom(frag, 1, Elements.C; name = "C\nD")  # newline inside atom name

    tmp = tempname() * ".cif"
    try
        write_mmcif(tmp, sys)
        text = read(tmp, String)
        atom_lines = filter(l -> startswith(l, "ATOM ") || startswith(l, "HETATM "), split(text, '\n'))
        @test length(atom_lines) == 1   # exactly one row, not split across lines
        # The newline must have been collapsed.
        @test !occursin('\n', atom_lines[1])
    finally
        rm(tmp; force=true)
    end
end
