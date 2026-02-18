@testitem "System" begin
    for T in [Float32, Float64]
        testsys = load_pdb(ball_data_path("../test/data/1tgh.pdb"), T)
        fdb = FragmentDB{T}()
        normalize_names!(testsys, fdb)
        reconstruct_fragments!(testsys, fdb)
        build_bonds!(testsys, fdb)

        @test natoms(testsys) == 3778
        @test nbonds(testsys) == 3839
        @test nmolecules(testsys) == 1
        @test nchains(testsys) == 6
        @test nsecondary_structures(testsys) == 30
        ct = chains(testsys)
        @test nfragments.(ct) == [12, 12, 180, 5, 1, 12]
        @test nfragments.(ct; variant = FragmentVariant.None) == [0, 0, 0, 5, 1, 12]
        @test nnucleotides.(ct) == [0, 0, 0, 0, 0, 0]
        @test nresidues.(ct) == [12, 12, 180, 0, 0, 0]

        # empty!
        sys = deepcopy(testsys)
        @test empty!(sys) === sys
        @test natoms(sys) == 0
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0
        @test nfragments(sys; variant = FragmentVariant.None) == 0
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 0

        # append!
        sys = deepcopy(testsys)
        @test append!(sys) === sys
        @test sys == testsys
        @test append!(sys, System{T}()) === sys
        @test sys == testsys
        @test append!(sys, testsys) === sys
        mol = first(molecules(sys))
        @test natoms(sys) == 2natoms(testsys)
        @test natoms(mol) == natoms(testsys)
        @test nbonds(sys) == 2nbonds(testsys)
        @test nbonds(mol) == nbonds(testsys)
        @test nmolecules(sys) == 2nmolecules(testsys)
        @test nchains(sys) == 2nchains(testsys)
        @test nchains(mol) == nchains(testsys)
        @test nsecondary_structures(sys) == 2nsecondary_structures(testsys)
        @test nsecondary_structures(mol) == nsecondary_structures(testsys)
        @test nfragments(sys) == 2nfragments(testsys)
        @test nfragments(mol) == nfragments(testsys)
        @test nfragments(sys; variant = FragmentVariant.None) == 2nfragments(testsys; variant = FragmentVariant.None)
        @test nfragments(mol; variant = FragmentVariant.None) == nfragments(testsys; variant = FragmentVariant.None)
        @test nnucleotides(sys) == 2nnucleotides(testsys)
        @test nnucleotides(mol) == nnucleotides(testsys)
        @test nresidues(sys) == 2nresidues(testsys)
        @test nresidues(mol) == nresidues(testsys)
        Atom(sys, 9999, Elements.N)
        @test allunique([
            atoms(sys).idx;
            bonds(sys).idx;
            molecules(sys).idx;
            chains(sys).idx;
            secondary_structures(sys).idx;
            fragments(sys).idx
        ])

        # atoms
        ## sort + sort! + sort_atoms!
        sys = deepcopy(testsys)
        at = atoms(sys)[100:2:200]
        at2 = sort(at; rev = true)
        @test at2 isa AtomTable{T}
        @test size(at2) == size(at)
        @test at2.idx == sort(at.idx; rev = true)

        at2 = sort(at; by = atom -> atom.name)
        @test at2 isa AtomTable{T}
        @test size(at2) == size(at)
        @test at2.name == sort(at.name)

        @test sort!(at; rev = true) === at
        @test at.idx == sort(atoms(sys)[100:2:200].idx; rev = true)

        @test sort!(at; by = atom -> atom.name) === at
        @test at.name == sort(atoms(sys)[100:2:200].name)

        @test sort_atoms!(sys; rev = true) === sys
        @test issorted(atoms(sys).idx; rev = true)

        @test sort_atoms!(sys; by = atom -> atom.name) === sys
        @test issorted(atoms(sys).name)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        at  = atoms(sys)[100:2:200]
        aidx = Set(at.idx)
        bidx = Set(bonds(at).idx)
        @test delete!(at) === at
        @test natoms(sys) == 3727
        @test length(at) == 0
        @test nbonds(sys) == 3721
        @test length(aidx ∩ Set(atoms(sys).idx)) == 0
        @test length(bidx ∩ Set(bonds(sys).idx)) == 0

        stale_table = atoms(sys)
        stale_col   = stale_table.idx

        at = atoms(sys)
        atom = first(at)
        aidx = atom.idx
        bidx = Set(bonds(atom).idx)
        @test_throws KeyError delete!(at, -1)
        @test delete!(at, aidx) === at
        @test natoms(sys) == 3726
        @test length(at) == 3726
        @test nbonds(sys) == 3719
        @test aidx ∉ atoms(sys).idx
        @test length(bidx ∩ Set(bonds(sys).idx)) == 0

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 3726
        @test length(stale_col) == 3726

        @test delete!(at) === at
        @test natoms(sys) == 0
        @test length(at) == 0
        @test nbonds(sys) == 0

        # bonds
        ## sort + sort! + sort_bonds!
        sys = deepcopy(testsys)
        bt = bonds(sys)[100:2:200]
        bt2 = sort(bt; rev = true)
        @test bt2 isa BondTable{T}
        @test size(bt2) == size(bt)
        @test bt2.idx == sort(bt.idx; rev = true)

        bt2 = sort(bt; by = bond -> bond.a2)
        @test bt2 isa BondTable{T}
        @test size(bt2) == size(bt)
        @test bt2.a2 == sort(bt.a2)

        @test sort!(bt; rev = true) === bt
        @test bt.idx == sort(bonds(sys)[100:2:200].idx; rev = true)

        @test sort!(bt; by = bond -> bond.a2) === bt
        @test bt.a2 == sort(bonds(sys)[100:2:200].a2)

        @test sort_bonds!(sys; rev = true) === sys
        @test issorted(bonds(sys).idx; rev = true)

        @test sort_bonds!(sys; by = bond -> bond.a2) === sys
        @test issorted(bonds(sys).a2)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        bt  = bonds(sys)[100:2:200]
        bidx = Set(bt.idx)
        @test delete!(bt) === bt
        @test length(bt) == 0
        @test nbonds(sys) == 3788
        @test length(bidx ∩ Set(bonds(sys).idx)) == 0

        stale_table = bonds(sys)
        stale_col   = stale_table.idx

        bt = bonds(sys)
        bidx = first(bt.idx)
        @test_throws KeyError delete!(bt, -1)
        @test delete!(bt, bidx) === bt
        @test length(bt) == 3787
        @test nbonds(sys) == 3787
        @test bidx ∉ bonds(sys).idx

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 3787
        @test length(stale_col) == 3787

        @test delete!(bt) === bt
        @test length(bt) == 0
        @test nbonds(sys) == 0

        # molecules
        ## sort + sort! + sort_molecules!
        sys = deepcopy(testsys)
        Molecule(sys)
        mt = molecules(sys)
        mt2 = sort(mt; rev = true)
        @test mt2 isa MoleculeTable{T}
        @test size(mt2) == size(mt)
        @test mt2.idx == sort(mt.idx; rev = true)

        mt2 = sort(mt; by = mol -> mol.name)
        @test mt2 isa MoleculeTable{T}
        @test size(mt2) == size(mt)
        @test mt2.name == sort(mt.name)

        @test sort!(mt; rev = true) === mt
        @test mt.idx == sort(molecules(sys).idx; rev = true)

        @test sort!(mt; by = mol -> mol.name) === mt
        @test mt.name == sort(molecules(sys).name)

        @test sort_molecules!(sys; rev = true) === sys
        @test issorted(molecules(sys).idx; rev = true)

        @test sort_molecules!(sys; by = mol -> mol.name) === sys
        @test issorted(molecules(sys).name)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        mt = molecules(sys)
        @test_throws KeyError delete!(mt, -1)
        @test delete!(mt, first(mt.idx); keep_atoms = true) === mt
        @test length(mt) == 0
        @test natoms(sys) == 3778
        @test all(isnothing, atoms(sys).molecule_idx)
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0

        sys = deepcopy(testsys)
        mt = molecules(sys)
        @test delete!(mt; keep_atoms = true) === mt
        @test length(mt) == 0
        @test natoms(sys) == 3778
        @test all(isnothing, atoms(sys).molecule_idx)
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0

        sys = deepcopy(testsys)
        mt = molecules(sys)
        stale_table = molecules(sys)
        stale_col   = stale_table.idx

        @test_throws KeyError delete!(mt, -1)
        @test delete!(mt, first(mt.idx); keep_atoms = false) === mt
        @test length(mt) == 0
        @test natoms(sys) == 0
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 0
        @test length(stale_col) == 0

        sys = deepcopy(testsys)
        mt = molecules(sys)
        @test delete!(mt; keep_atoms = false) === mt
        @test length(mt) == 0
        @test natoms(sys) == 0
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0

        # chains
        ## sort + sort! + sort_chains!
        sys = deepcopy(testsys)
        ct = chains(sys)[1:2]
        ct2 = sort(ct; rev = true)
        @test ct2 isa ChainTable{T}
        @test size(ct2) == size(ct)
        @test ct2.idx == sort(ct.idx; rev = true)

        ct2 = sort(ct; by = chain -> chain.name)
        @test ct2 isa ChainTable{T}
        @test size(ct2) == size(ct)
        @test ct2.name == sort(ct.name)

        @test sort!(ct; rev = true) === ct
        @test ct.idx == sort(chains(sys)[1:2].idx; rev = true)

        @test sort!(ct; by = chain -> chain.name) === ct
        @test ct.name == sort(chains(sys)[1:2].name)

        @test sort_chains!(sys; rev = true) === sys
        @test issorted(chains(sys).idx; rev = true)

        @test sort_chains!(sys; by = chain -> chain.name) === sys
        @test issorted(chains(sys).name)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        ct = chains(sys)[1:2]
        @test delete!(ct; keep_atoms = true) === ct
        @test length(ct) == 0
        @test natoms(sys) == 3778
        @test natoms(sys; chain_idx = Some(nothing)) == 768
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 4
        @test nsecondary_structures(sys) == 28
        @test nfragments(sys) == 198

        ct = chains(sys)
        @test delete!(ct; keep_atoms = false) === ct
        @test length(ct) == 0
        @test natoms(sys) == 768
        @test nbonds(sys) == 822
        @test nmolecules(sys) == 1
        @test nchains(sys) == 0
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 0

        sys = deepcopy(testsys)
        stale_table = chains(sys)
        stale_col   = stale_table.idx

        ct = chains(sys)
        @test_throws KeyError delete!(ct, -1)
        @test delete!(ct, first(ct.idx); keep_atoms = true) === ct
        @test length(ct) == 5
        @test natoms(sys) == 3778
        @test natoms(sys; chain_idx = Some(nothing)) == 384
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 5
        @test nsecondary_structures(sys) == 29
        @test nfragments(sys) == 210

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 5
        @test length(stale_col) == 5

        @test delete!(ct, first(ct.idx); keep_atoms = false) === ct
        @test length(ct) == 4
        @test natoms(sys) == 3394
        @test nbonds(sys) == 3428
        @test nmolecules(sys) == 1
        @test nchains(sys) == 4
        @test nsecondary_structures(sys) == 28
        @test nfragments(sys) == 198

        # secondary structures
        ## sort + sort! + sort_secondary_structures!
        sys = deepcopy(testsys)
        st = secondary_structures(sys)[1:2]
        st2 = sort(st; rev = true)
        @test st2 isa SecondaryStructureTable{T}
        @test size(st2) == size(st)
        @test st2.idx == sort(st.idx; rev = true)

        st2 = sort(st; by = ss -> ss.name)
        @test st2 isa SecondaryStructureTable{T}
        @test size(st2) == size(st)
        @test st2.name == sort(st.name)

        @test sort!(st; rev = true) === st
        @test st.idx == sort(secondary_structures(sys)[1:2].idx; rev = true)

        @test sort!(st; by = ss -> ss.name) === st
        @test st.name == sort(secondary_structures(sys)[1:2].name)

        @test sort_secondary_structures!(sys; rev = true) === sys
        @test issorted(secondary_structures(sys).idx; rev = true)

        @test sort_secondary_structures!(sys; by = ss -> ss.name) === sys
        @test issorted(secondary_structures(sys).name)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        st = secondary_structures(sys)[1:2]
        @test delete!(st; keep_fragments = true) === st
        @test length(st) == 0
        @test natoms(sys) == 3778
        @test natoms(sys) - natoms(secondary_structures(sys)) == 822
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 28
        @test nfragments(sys) == 222

        st = secondary_structures(sys)
        @test delete!(st; keep_fragments = false) === st
        @test length(st) == 0
        @test natoms(sys) == 822
        @test nbonds(sys) == 858
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 0
        @test nfragments(sys) == 42

        sys = deepcopy(testsys)
        stale_table = secondary_structures(sys)
        stale_col   = stale_table.idx

        st = secondary_structures(sys)
        @test_throws KeyError delete!(st, -1)
        @test delete!(st, first(st.idx); keep_fragments = true) === st
        @test length(st) == 29
        @test natoms(sys) == 3778
        @test natoms(sys) - natoms(secondary_structures(sys)) == 438
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 29
        @test nfragments(sys) == 222

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 29
        @test length(stale_col) == 29

        @test delete!(st, first(st.idx); keep_fragments = false) === st
        @test length(st) == 28
        @test natoms(sys) == 3394
        @test nbonds(sys) == 3428
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 28
        @test nfragments(sys) == 210

        # fragments
        ## sort + sort! + sort_fragments!
        sys = deepcopy(testsys)
        ft = fragments(sys)[1:2:end]
        ft2 = sort(ft; rev = true)
        @test ft2 isa FragmentTable{T}
        @test size(ft2) == size(ft)
        @test ft2.idx == sort(ft.idx; rev = true)

        ft2 = sort(ft; by = frag -> frag.name)
        @test ft2 isa FragmentTable{T}
        @test size(ft2) == size(ft)
        @test ft2.name == sort(ft.name)

        @test sort!(ft; rev = true) === ft
        @test ft.idx == sort(fragments(sys)[1:2:end].idx; rev = true)

        @test sort!(ft; by = frag -> frag.name) === ft
        @test ft.name == sort(fragments(sys)[1:2:end].name)

        @test sort_fragments!(sys; rev = true) === sys
        @test issorted(fragments(sys).idx; rev = true)

        @test sort_fragments!(sys; by = frag -> frag.name) === sys
        @test issorted(fragments(sys).name)

        ## delete! + revalidate_indices!
        sys = deepcopy(testsys)
        ft = fragments(sys)[1:2:end]
        @test delete!(ft; keep_atoms = true) === ft
        @test length(ft) == 0
        @test natoms(sys) == 3778
        @test natoms(sys; fragment_idx = Some(nothing)) == 1917
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 30
        @test nfragments(sys) == 111

        ft = fragments(sys)
        @test delete!(ft; keep_atoms = false) === ft
        @test length(ft) == 0
        @test natoms(sys) == 1917
        @test nbonds(sys) == 1843
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 30
        @test nfragments(sys) == 0

        sys = deepcopy(testsys)
        stale_table = fragments(sys)
        stale_col   = stale_table.idx

        ft = fragments(sys)
        @test_throws KeyError delete!(ft, -1)
        @test delete!(ft, first(ft.idx); keep_atoms = true) === ft
        @test length(ft) == 221
        @test natoms(sys) == 3778
        @test natoms(sys; fragment_idx = Some(nothing)) == 31
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 30
        @test nfragments(sys) == 221

        @test_throws KeyError first(stale_table.idx)
        @test_throws KeyError first(stale_col)
        @test revalidate_indices!(stale_table) === stale_table
        @test revalidate_indices!(stale_col) === stale_col
        @test length(stale_table) == 221
        @test length(stale_col) == 221

        @test delete!(ft, first(ft.idx); keep_atoms = false) === ft
        @test length(ft) == 220
        @test natoms(sys) == 3745
        @test nbonds(sys) == 3802
        @test nmolecules(sys) == 1
        @test nchains(sys) == 6
        @test nsecondary_structures(sys) == 30
        @test nfragments(sys) == 220
    end
end
