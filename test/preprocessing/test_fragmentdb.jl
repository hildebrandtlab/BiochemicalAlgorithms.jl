@testitem "FragmentDB" begin
    default_fdb = FragmentDB()

    @test length(default_fdb.fragments) == 33
    @test length(default_fdb.name_mappings) == 6
    @test length(default_fdb.defaults) == 1

    fdb = FragmentDB(ball_data_path("fragments/Fragments.db.json"))

    @test fdb == default_fdb
end

@testitem "normalize_names!()" begin
    default_fdb = FragmentDB()

    # this file contains some manually changed names, e.g. "ATOM     13 1HB  ARG ..." is changed to "ATOM     13 HB2  ARG ..."
    sys = load_pdb(ball_data_path("../test/data/AmberFF_naming.pdb"))
    @test natoms(sys) == 892

    # check Amber naming 
    arg = fragments(sys)[1]
    @test atoms(arg)[13].name == "HB2"
    @test atoms(arg)[19].name == "HB3"

    # check Charmm naming
    @test atoms(arg)[16].name == "HH1"
    @test atoms(arg)[17].name == "HH2"

    
    # check amber
    @test atoms(arg)[21].name == "HG2"

    normalize_names!(sys, FragmentDB())


    @test atoms(arg)[13].name == "1HB"
    @test atoms(arg)[19].name == "2HB"

    # check Charmm naming
    @test atoms(arg)[16].name == "1HH1"
    @test atoms(arg)[17].name == "1HH2"

    # check amber (if it would be discover than "2HG")
    @test atoms(arg)[21].name == "1HG"
end

#TODO: Should we check all the naming conventions
#TODO: Add test cases for reconstruct_fragments!() and 
#TODO: Add test cases for build_bonds!()
