@testset "FragmentDB" begin
    default_fdb = FragmentDB()

    @test length(default_fdb.fragments) == 33
    @test length(default_fdb.name_mappings) == 6
    @test length(default_fdb.defaults) == 1

    fdb = FragmentDB(ball_data_path("fragments/Fragments.db.json"))

    @test fdb == default_fdb
end