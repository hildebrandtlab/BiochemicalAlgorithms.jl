@testitem "SystemComponentTableCol" begin
    using BiochemicalAlgorithms: SystemComponentTableCol

    for T in [Int, Float32, Float64]
        v = T[1, 2, 3]

        # invalid row bounds
        @test_throws BoundsError SystemComponentTableCol(v, [0], Dict(0 => 0))[1]
        @test_throws BoundsError SystemComponentTableCol(v, [4], Dict(4 => 4))[1]

        # invalid mapping
        @test_throws KeyError SystemComponentTableCol(v, [1], Dict(0 => 0))[1]

        # row subset
        u = SystemComponentTableCol(v, [3, 1], Dict(3 => 3, 1 => 1))
        @test u isa SystemComponentTableCol{T}
        @test eltype(u) === T
        @test size(u) == (2,)
        @test length(u) == 2
        @test_throws BoundsError u[0]
        @test_throws BoundsError u[0] = T(41)
        @test u[1] == v[3]
        u[1] = 11
        @test u[1] == v[3] == 11
        @test u[2] == v[1]
        @test_throws BoundsError u[3]
        @test_throws BoundsError u[3] = T(42)

        u = SystemComponentTableCol(v, [3, 1], Dict(3 => 1, 1 => 3))
        @test u[1] == v[1]
        u[1] = 11
        @test u[1] == v[1] == 11
        @test u[2] == v[3]

        # row bag
        u = SystemComponentTableCol(v, [2, 2], Dict(2 => 2))
        @test size(u) == (2,)
        @test length(u) == 2
        @test_throws BoundsError u[0]
        @test_throws BoundsError u[0] = T(41)
        @test u[1] == u[2] == v[2]
        u[1] = 11
        @test u[1] == u[2] == v[2] == 11
        u[2] = 12
        @test u[1] == u[2] == v[2] == 12
        @test_throws BoundsError u[3]
        @test_throws BoundsError u[3] = T(42)
    end
end
