@testitem "FragmentDB" begin
    for fdb in (FragmentDB(), FragmentDB{Float32}(), FragmentDB{Float64}())
        @test length(fdb.fragments) == 33
        @test length(fdb.name_mappings) == 6
        @test length(fdb.defaults) == 1
    end
end
