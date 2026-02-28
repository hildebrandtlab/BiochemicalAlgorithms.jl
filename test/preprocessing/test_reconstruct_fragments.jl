@testitem "reconstruct_fragments!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # single residue (LYS)
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="LYS")

            a = Atom(f, 1, Elements.C; name="CA")

            num_reconstructed = reconstruct_fragments!(sys, fdb)

            @test num_reconstructed == natoms(sys) - 1
            @test natoms(f) == length(fdb.fragments["LYS"].variants[1].atoms)
            @test natoms(f) == 22
        end

        let sys = load_hinfile(ball_data_path("../test/data/ReconstructFragmentProcessor_test1.hin"), T)
            normalize_names!(sys, fdb)

            @test reconstruct_fragments!(sys, fdb) == 4
            @test natoms(sys) == 31
        end

        # all 20 standard amino acids: single CA atom → reconstruct → verify atom count
        expected_atoms = Dict(
            "ALA" => 10, "ARG" => 24, "ASN" => 14, "ASP" => 12, "CYS" => 11,
            "GLN" => 17, "GLU" => 15, "GLY" =>  7, "HIS" => 18, "ILE" => 19,
            "LEU" => 19, "LYS" => 22, "MET" => 17, "PHE" => 20, "PRO" => 14,
            "SER" => 11, "THR" => 14, "TRP" => 24, "TYR" => 21, "VAL" => 16,
        )
        for (aa, expected) in expected_atoms
            sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name=aa)
            Atom(f, 1, Elements.C; name="CA")

            n = reconstruct_fragments!(sys, fdb)
            @test n == expected - 1
            @test natoms(f) == expected
        end

        # multi-fragment system: already complete → 0 reconstructed
        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            normalize_names!(sys, fdb)
            before = natoms(sys)
            @test reconstruct_fragments!(sys, fdb) == 0
            @test natoms(sys) == before
        end

        # unknown fragment: warning, 0 atoms added, no crash
        let sys = System{T}()
            f = Fragment(Chain(Molecule(sys)), 1; name="UNKNOWN")
            Atom(f, 1, Elements.C; name="CA")
            @test reconstruct_fragments!(sys, fdb) == 0
            @test natoms(f) == 1
        end

        # empty system: 0 reconstructed
        let sys = System{T}()
            @test reconstruct_fragments!(sys, fdb) == 0
        end
    end
end
