@testitem "Bounding box" begin
    for T in [Float32, Float64]
        # empty system
        let sys = System{T}()
            box = compute_bounding_box(sys)
            @test box isa BoundingBox{T}
            @test box.min == zero(Vector3{T})
            @test box.max == zero(Vector3{T})
        end

        let sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
            box = compute_bounding_box(sys)
            @test box isa BoundingBox{T}
            @test box.min ≈ T[-1.801, -1.358, -1.827]
            @test box.max ≈ T[1.676, 3.039, 5.473]

            box = compute_bounding_box(first(fragments(sys)))
            @test box isa BoundingBox{T}
            @test box.min ≈ T[-1.801, -0.928, -1.827]
            @test box.max ≈ T[1.642, 2.118, 1.838]
        end
    end
end
