@testitem "add_hydrogens!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # AlaAla dipeptide: adds terminal OXT hydrogen
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            normalize_names!(sys, fdb)
            reconstruct_fragments!(sys, fdb)
            build_bonds!(sys, fdb)

            before = natoms(sys)
            @test add_hydrogens!(sys) == 1
            @test natoms(sys) == before + 1

            # verify the added atom is a hydrogen
            added = last(atoms(sys))
            @test added.element == Elements.H
            @test added.name == "HOXT"
        end

        # AlaGlySer tripeptide: adds terminal OXT hydrogen
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)
            reconstruct_fragments!(sys, fdb)
            build_bonds!(sys, fdb)

            @test add_hydrogens!(sys) == 1
            @test last(atoms(sys)).name == "HOXT"
        end

        # bpti: adds terminal OXT hydrogen
        let sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
            normalize_names!(sys, fdb)
            reconstruct_fragments!(sys, fdb)
            build_bonds!(sys, fdb)

            @test add_hydrogens!(sys) == 1
            @test last(atoms(sys)).name == "HOXT"
        end

        # already-saturated system: second pass adds 0
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            normalize_names!(sys, fdb)
            reconstruct_fragments!(sys, fdb)
            build_bonds!(sys, fdb)

            add_hydrogens!(sys)
            n = natoms(sys)
            @test add_hydrogens!(sys) == 0
            @test natoms(sys) == n
        end

        # empty system: no crash
        let sys = System{T}()
            @test add_hydrogens!(sys) == 0
        end
    end
end
