@testitem "Vector3" begin
    for T in [Int, Float32, Float64]
        u = Vector3{T}([-1, -2, -3])
        v = Vector3{T}(1, 2, 3)

        @test v isa Vector3{T}
        @test eltype(v) == T
        @test size(v) == (3,)
        @test length(v) == 3
        @test v == T[1, 2, 3]
        @test v[1] == 1
        @test v[2] == 2
        @test v[3] == 3
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[4]
        @test_throws DimensionMismatch Vector3{T}()

        @test zeros(Vector3{T}) isa Vector3{T}
        @test zeros(Vector3{T}) == T[0, 0, 0]
        @test ones(Vector3{T}) isa Vector3{T}
        @test ones(Vector3{T}) == T[1, 1, 1]

        @test squared_norm(v) ≈ 1 * 1 + 2 * 2 + 3 * 3
        @test distance(u, v) == distance(v, u) == 2sqrt(squared_norm(v))
    end
end

@testitem "Matrix3" begin
    for T in [Int, Float32, Float64]
        m = Matrix3{T}((1, 2, 3, 4, 5, 6, 7, 8, 9))
        @test m isa Matrix3{T}
        @test eltype(m) == T
        @test size(m) == (3, 3)
        @test length(m) == 9
        t = T[1 4 7; 2 5 8; 3 6 9]
        @test m == t
        j = 1
        for i in m
            @test i == t[j]
            j += 1
        end
        @test_throws BoundsError m[0]
        @test_throws BoundsError m[10]
        @test_throws DimensionMismatch Matrix3{T}()
    end
end

@testitem "Properties" begin
    @test Properties() isa Properties
    @test Properties(:A => 23, :B => 12, :C => 12) isa Properties
end

@testitem "Flags" begin
    @test Flags() isa Flags
    @test Flags([:first, :second]) isa Flags
end

@testitem "RowProjectionVector" begin
    using BiochemicalAlgorithms: RowProjectionVector

    for T in [Int, Float32, Float64]
        v = T[1, 2, 3]
        
        # invalid row bounds
        @test_throws BoundsError RowProjectionVector(v, [0])[1]
        @test_throws BoundsError RowProjectionVector(v, [4])[1]

        # row subset
        u = RowProjectionVector(v, [3, 1])
        @test u isa RowProjectionVector{T}
        @test size(u) == (2,)
        @test length(u) == 2
        @test_throws BoundsError u[0]
        @test_throws BoundsError u[0] = T(41)
        @test u[1] == v[3]
        u[1] = T(11)
        @test u[1] == v[3] == 11
        @test u[2] == v[1]
        @test_throws BoundsError u[3]
        @test_throws BoundsError u[3] = T(42)

        # row bag
        u = RowProjectionVector(v, [2, 2])
        @test size(u) == (2,)
        @test length(u) == 2
        @test_throws BoundsError u[0]
        @test_throws BoundsError u[0] = T(41)
        @test u[1] == u[2] == v[2]
        u[1] = T(11)
        @test u[1] == u[2] == v[2] == 11
        u[2] = T(12)
        @test u[1] == u[2] == v[2] == 12
        @test_throws BoundsError u[3]
        @test_throws BoundsError u[3] = T(42)
    end
end
