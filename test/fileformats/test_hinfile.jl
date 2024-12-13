@testitem "reading HIN files" begin
    for T in [Float32, Float64]
        sys = load_hinfile(ball_data_path("../test/data/hinfile_test.hin"), T)

        @test natoms(sys) == 648
        @test nmolecules(sys) == 216
        @test nfragments(sys) == 0

        a = first(atoms(sys))
        @test a.name == "O"
        @test a.element == Elements.O
        @test a.charge ≈ -0.834
        @test a.r ≈ Vector3{T}(0.5903809, -0.4102751, -0.8605154)
        @test nbonds(a) == 2

        @test get_property(sys, :temperature) ≈ 297.5626

        @test get_property(sys, :periodic_box_width)  ≈ 18.70136
        @test get_property(sys, :periodic_box_height) ≈ 18.70136
        @test get_property(sys, :periodic_box_depth)  ≈ 18.70136

        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)

        @test natoms(sys) == 31
        @test nmolecules(sys) == 1
        @test nfragments(sys) == 3
        @test nresidues(sys) == 3
        @test nchains(sys) == 1
        @test nbonds(sys) == 30

        @test_throws SystemError load_hinfile(ball_data_path("../test/data/ASDFASDFASEFADSFASDFAEW.hin"), T)
        @test_throws MethodError load_hinfile(ball_data_path("../test/data/hinfile_test_invalid.hin"), T)
    end
end

@testitem "writing HIN files" begin
    for T in [Float32, Float64]
        sys = load_hinfile(ball_data_path("../test/data/hinfile_test.hin"), T)

        outfname = tempname()
        write_hinfile(outfname, sys)

        sys2 = load_hinfile(outfname, T)

        @test sys == sys2

        # test name truncation
        first(atoms(sys)).name = "TEST NAME"

        write_hinfile(outfname, sys)

        sys2 = load_hinfile(outfname, T)

        @test first(atoms(sys2)).name == "TEST"
    end
end
