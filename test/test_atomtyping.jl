@testset "DEF_file_parser" begin
    def_file = load_atomtyping_DEF("data/TEST_DEF_file_parser_ATOMTYPE_GFF.DEF")

    @test def_file isa DataFrame
    @test nrow(def_file) == 296
    @test ncol(def_file) == 8
    @test def_file.type_name[9] == "ca"
    @test all(in(["*"]).(def_file.residue_names))
    @test typeof(def_file.atomic_number) == Vector{Int64}
    @test def_file.atomic_property[5] == "[1DB,0DL]"
    @test def_file.CES[18] == "(XD3[sb',db])"

end
