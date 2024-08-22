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
