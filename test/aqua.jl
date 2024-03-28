@testitem "Aqua" begin
    using Aqua

    # Current workaround to avoid failing tests due to ambiguity in dependencies
    # https://github.com/JuliaTesting/Aqua.jl/issues/77
    @testset "Method ambiguity" begin
        Aqua.test_ambiguities(BiochemicalAlgorithms)
    end

    Aqua.test_all(BiochemicalAlgorithms;
        ambiguities=false
    )
end
