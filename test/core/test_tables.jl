@testitem "AbstractColumnTable" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        at  = atoms(sys)
        buf = IOBuffer()

        show(buf, at)
        @test !isempty(take!(buf))

        show(buf, MIME"text/html"(), at)
        @test !isempty(take!(buf))

        show(buf, MIME"text/markdown"(), at)
        @test !isempty(take!(buf))

        show(buf, MIME"text/plain"(), at)
        @test !isempty(take!(buf))
    end
end
