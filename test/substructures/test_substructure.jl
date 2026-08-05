@testitem "Substructure" begin
    for T in [Float32, Float64]
        # read a simple molecule
        sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)
        filter_fn = atom -> atom.element == Elements.C

        mol = first(molecules(load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)))

        # create a substructure
        filtered_atoms = filter(filter_fn, atoms(mol))
        idxset = Set(filtered_atoms.idx)
        filtered_bonds = filter(bond ->
            bond.atom1_idx ∈ idxset && bond.atom2_idx ∈ idxset,
            bonds(mol)
        )

        s = Substructure(
            "aspirin substructure",
            mol,
            filtered_atoms,
            filtered_bonds
        )

        @test s isa Substructure{T}
        @test s.name == "aspirin substructure"
        @test natoms(s) == 9
        @test nbonds(s) == 8
        @test length(s.properties) == 3

        # generate the same substructure by filtering
        s2 = filter_atoms(
            filter_fn,
            mol;
            name="aspirin substructure",
            adjacent_bonds=false
        )

        @test s2 == s

        sys = copy(s)
        @test sys.name == "aspirin substructure"
        @test length(sys.properties) == 3
        @test natoms(sys) == natoms(s)
        @test nbonds(sys) == nbonds(s)
        @test nmolecules(sys) == nmolecules(s)
        @test nproteins(sys) == nproteins(s)
        @test nchains(sys) == nchains(s)
        @test nfragments(sys) == nfragments(s)
        @test nresidues(sys) == nresidues(s)
        @test nnucleotides(sys) == nnucleotides(s)
        @test nsecondary_structures(sys) == nsecondary_structures(s)

        # empty substructure (filter matches nothing)
        s_empty = filter_atoms(a -> false, mol; name="empty", adjacent_bonds=false)
        @test natoms(s_empty) == 0
        @test nbonds(s_empty) == 0

        sys_empty = copy(s_empty)
        @test sys_empty isa System{T}
        @test natoms(sys_empty) == 0
        @test nbonds(sys_empty) == 0
        @test nmolecules(sys_empty) == 0

        # check substruct interface for all atom containers
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        fdb = FragmentDB{T}()
        infer_topology!(sys, fdb)

        for ac in (sys, first(molecules(sys)), first(chains(sys)), first(fragments(sys)), first(bonds(sys)))
            filtered_atoms = filter(a -> true, atoms(ac))
            idxset = Set(filtered_atoms.idx)
            filtered_bonds = filter(bond ->
                bond.atom1_idx ∈ idxset && bond.atom2_idx ∈ idxset,
                bonds(parent(ac))
            )

            s = Substructure(
                "AlaAla",
                ac,
                filtered_atoms,
                filtered_bonds
            )

            @test s isa Substructure{T}
            @test s.name == "AlaAla"
            @test natoms(s) > 0
            @test nbonds(s) > 0
            @test nmolecules(s) > 0
            @test nproteins(s) == 0
            @test nchains(s) > 0
            @test nfragments(s) > 0
            @test nresidues(s) > 0
            @test nnucleotides(s) == 0
            @test nsecondary_structures(s) > 0

            # generate the same substructure by filtering
            s2 = filter_atoms(
                a -> true,
                ac;
                name="AlaAla",
                adjacent_bonds=false
            )

            @test s2 == s
        end
    end
end
