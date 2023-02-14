@testset "Fragment" begin
    # simple fragment
    f = ( number = 1, name = "my fragment", chain = "protein chain A")
    @test f isa Fragment
    @test f.number == 1
    @test f.name == "my fragment"
    @test f.chain == "protein chain A"

    # incomplete fragment
    f2 = ( number = 1, name = "my fragment")
    @test !isa(f2, Fragment)
end
