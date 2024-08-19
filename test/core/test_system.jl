@testitem "System" begin
    for T in [Float32, Float64]
        testsys = load_pdb(ball_data_path("../test/data/1tgh.pdb"), T)
        fdb = FragmentDB{T}()
        normalize_names!(testsys, fdb)
        reconstruct_fragments!(testsys, fdb)
        build_bonds!(testsys, fdb)

        @test natoms(testsys) == 3790
        @test nbonds(testsys) == 3839
        @test nmolecules(testsys) == 1
        @test nchains(testsys) == 3
        ct = chains(testsys)
        @test nfragments.(ct) == [192, 17, 13]
        @test nfragments.(ct; variant = FragmentVariant.None) == [12, 5, 1]
        @test nnucleotides.(ct) == [0, 12, 12]
        @test nresidues.(ct) == [180, 0, 0]

        # empty!
        sys = deepcopy(testsys)
        @test empty!(sys) === sys
        @test natoms(sys) == 0
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nfragments(sys) == 0
        @test nfragments(sys; variant = FragmentVariant.None) == 0
        @test nnucleotides(sys) == 0
        @test nresidues(sys) == 0

        # delete! atoms
        sys = deepcopy(testsys)
        at  = atoms(sys)[100:2:200]
        aidx = Set(at.idx)
        bidx = Set(bonds(at).idx)
        @test delete!(at) === at
        @test natoms(sys) == 3790 - length(aidx)
        @test length(at) == 0
        @test nbonds(sys) == 3839 - length(bidx)
        @test length(aidx ∩ Set(atoms(sys).idx)) == 0
        @test length(bidx ∩ Set(bonds(sys).idx)) == 0

        at = atoms(sys)
        delete!(at) === at
        @test natoms(sys) == 0
        @test length(at) == 0
        @test nbonds(sys) == 0

        # delete! bonds
        sys = deepcopy(testsys)
        bt  = bonds(sys)[100:2:200]
        bidx = Set(bt.idx)
        @test delete!(bt) === bt
        @test length(bt) == 0
        @test nbonds(sys) == 3839 - length(bidx)
        @test length(bidx ∩ Set(bonds(sys).idx)) == 0

        bt = bonds(sys)
        delete!(bt) === bt
        @test length(bt) == 0
        @test nbonds(sys) == 0

        # delete! molecules
        sys = deepcopy(testsys)
        mt = molecules(sys)
        delete!(mt; keep_atoms = true) === mt
        @test length(mt) == 0
        @test natoms(sys) == 3790
        @test all(isnothing, atoms(sys).molecule_idx)
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nfragments(sys) == 0

        sys = deepcopy(testsys)
        mt = molecules(sys)
        delete!(mt; keep_atoms = false) === mt
        @test length(mt) == 0
        @test natoms(sys) == 0
        @test nbonds(sys) == 0
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nfragments(sys) == 0

        # delete! chains
        sys = deepcopy(testsys)
        ct = chains(sys)[1:2]
        delete!(ct; keep_atoms = true) === ct
        @test length(ct) == 0
        @test natoms(sys) == 3790
        @test natoms(sys; chain_idx = Some(nothing)) == 3401
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 13

        ct = chains(sys)
        delete!(ct; keep_atoms = false) === ct
        @test length(ct) == 0
        @test natoms(sys) == 3401
        @test nbonds(sys) == 3839 - 413
        @test nmolecules(sys) == 1
        @test nchains(sys) == 0
        @test nfragments(sys) == 0

        # delete! fragments
        sys = deepcopy(testsys)
        ft = fragments(sys)[1:2:end]
        delete!(ft; keep_atoms = true) === ft
        @test length(ft) == 0
        @test natoms(sys) == 3790
        @test natoms(sys; fragment_idx = Some(nothing)) == 1929
        @test nbonds(sys) == 3839
        @test nmolecules(sys) == 1
        @test nchains(sys) == 3
        @test nfragments(sys) == 111

        ft = fragments(sys)
        delete!(ft; keep_atoms = false) === ft
        @test length(ft) == 0
        @test natoms(sys) == 1929
        @test nbonds(sys) == 3839 - 1984
        @test nmolecules(sys) == 1
        @test nchains(sys) == 3
        @test nfragments(sys) == 0
    end
end
