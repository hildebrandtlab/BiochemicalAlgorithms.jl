@testitem "build_bonds!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            empty!(sys._bonds)

            # normalize_names! is still required for terminal labeling
            infer_topology!(sys, fdb; reconstruct_fragments = false)
            @test nbonds.(fragments(sys)) == [12, 8, 12]
        end

        let sys = load_pdb(ball_data_path("../test/data/AmberFF_bench.pdb"), T)
            empty!(sys._bonds)

            infer_topology!(sys, fdb; reconstruct_fragments = false)
            @test nbonds(sys) == 906
        end

        let sys = load_pdb(ball_data_path("../test/data/ACE_test_A.pdb"), T)
            infer_topology!(sys, fdb; reconstruct_fragments = false)
            @test nbonds(sys) == 1666
        end

        let sys = load_pdb(ball_data_path("../test/data/ACE_test_B.pdb"), T)
            infer_topology!(sys, fdb; reconstruct_fragments = false)
            @test nbonds(sys) == 468

            # verify disulfide bridges (BPTI has 3: CYS5-CYS55, CYS14-CYS38, CYS30-CYS51)
            ss_bonds = filter(b -> :TYPE__DISULPHIDE_BOND in b.flags, bonds(sys))
            @test length(ss_bonds) == 3

            for b in ss_bonds
                a1 = atom_by_idx(sys, b.a1)
                a2 = atom_by_idx(sys, b.a2)
                @test a1.name == "SG"
                @test a2.name == "SG"
                @test parent_fragment(a1).name == "CYS"
                @test parent_fragment(a2).name == "CYS"
            end
        end

        # bond order verification on AlaAla dipeptide
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            infer_topology!(sys, fdb)

            @test count(b -> b.order == BondOrder.Single, bonds(sys)) == 20
            @test count(b -> b.order == BondOrder.Double, bonds(sys)) == 2

            # verify inter-fragment peptide bond
            peptide_bonds = filter(bonds(sys)) do b
                a1 = atom_by_idx(sys, b.a1)
                a2 = atom_by_idx(sys, b.a2)
                parent_fragment(a1) != parent_fragment(a2)
            end
            @test length(peptide_bonds) == 1
            @test peptide_bonds[1].order == BondOrder.Single
        end

        # idempotency: building bonds twice should not create duplicates
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            infer_topology!(sys, fdb)
            n = nbonds(sys)
            @test n > 0

            build_bonds!(sys, fdb)
            @test nbonds(sys) == n
        end

        # empty system
        let sys = System{T}()
            build_bonds!(sys, fdb)
            @test nbonds(sys) == 0
        end

        # unknown fragment
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="ZZZZZ")
            Atom(f, 1, Elements.C; name="CA")
            build_bonds!(sys, fdb)
            @test nbonds(sys) == 0
        end
    end
end
