@testitem "reconstruct_fragments!()" begin
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        # single residue
        s = System{T}()
        f = Fragment(Chain(Molecule(s)), 1; name="LYS")

        a = Atom(f, 1, Elements.C; name="CA")

        num_reconstructed = reconstruct_fragments!(s, fdb)

        @test num_reconstructed == natoms(s) - 1
        @test natoms(f) == length(fdb.fragments["LYS"].variants[1].atoms)
        @test natoms(f) == 22
    end
end