@testset "Fragment" begin
    # simple fragment
    f = (idx = 0, number = 1, name = "my fragment", properties = Properties(), flags = Flags())

    @test f isa FragmentTuple
    @test f.number == 1
    @test f.name == "my fragment"
end
