@testitem "Aqua" begin
    using Aqua

    # Ignore hotfixed package IntervalSets (https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/pull/221)
    # for downgrade checks where the package is not pulled in (e.g. via MakieCore)
    Aqua.test_all(BiochemicalAlgorithms; stale_deps = (; ignore = [:IntervalSets]))
end
