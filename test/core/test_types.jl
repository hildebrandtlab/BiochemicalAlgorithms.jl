@testitem "Simple Vector3" begin
    v = Vector3{Int}(1, 2, 3)

    @test v isa Vector3{Int}
    @test size(v,1) == 3
    @test v[1] == 1
    @test v[2] == 2
    @test v[3] == 3
    @test_throws BoundsError v[4]
    @test_throws DimensionMismatch Vector3{Int}()
end


@testitem "Matrix3" begin
    m = Matrix3{Int}((1, 2, 3, 4, 5, 6, 7, 8, 9))
    @test m isa Matrix3{Int}
    @test size(m) == (3, 3)
    t = [1 4 7; 2 5 8; 3 6 9]
    j = 1
    for i in m
        @test i == t[j]
        global j += 1
    end
    @test_throws BoundsError m[10]
    @test_throws DimensionMismatch Matrix3{Int}()
end

@testitem "Properties" begin
    p = Properties()
    @test p isa Properties

    p2 = Properties(:A => 23, :B => 12, :C => 12)
    @test p2 isa Properties
end
