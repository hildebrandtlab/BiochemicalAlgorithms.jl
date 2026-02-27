@testitem "reconstruct_fragments!" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # single residue
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
    end
end
