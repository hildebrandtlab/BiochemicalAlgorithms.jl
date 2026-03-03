@testitem "Constants" begin
    using BiochemicalAlgorithms.Constants

    # vacuum electric permittivity (CODATA 2018)
    @test ε₀_nounit isa Float64
    @test isapprox(ε₀_nounit, 8.8541878128e-12; rtol=1e-10)

    # Avogadro's number (CODATA 2018)
    @test N_A_nounit isa Float64
    @test isapprox(N_A_nounit, 6.02214076e23; rtol=1e-10)

    # elementary charge (CODATA 2018)
    @test e₀_nounit isa Float64
    @test isapprox(e₀_nounit, 1.602176634e-19; rtol=1e-10)
end
