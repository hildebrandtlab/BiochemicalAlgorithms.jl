function _compare_without_system(m1::AbstractMolecule, m2::AbstractMolecule)
    result =       m1.name == m2.name &&
             atoms_df(m1)  == atoms_df(m2) &&
             bonds_df(m1)  == bonds_df(m2) &&
             m1.properties == m2.properties

    result
end

@testset "SDFile" begin
    @testset "Reading" begin
        mols = molecules(load_sdfile("data/sdfile_test_1.sdf"))

        @test length(mols) == 11

        @test sum(map(natoms, mols)) == 518
        @test sum(map(nbonds, mols)) == 528

        mol = mols[1]
        @test natoms(mol) == 39
        @test nbonds(mol) == 42

    @test has_property(mol, :NAME)
    @test get_property(mol, :NAME) == "Abacavir_sulfate"

    @test has_property(mol, :b_1rotN)
    @test get_property(mol, :b_1rotN) == "6"

    @test has_property(mol, :Weight)
    @test get_property(mol, :Weight) == "286.339"

    @test has_property(mol, :TPSA)
    @test get_property(mol, :TPSA) == "101.88"

    @test has_property(mol, :a_acc)
    @test get_property(mol, :a_acc) == "4"

    @test has_property(mol, :a_don)
    @test get_property(mol, :a_don) == "3"

    @test has_property(mol, Symbol("logP(o/w)"))
    @test get_property(mol, Symbol("logP(o/w)")) == "0.40906"

    @test has_property(mol, :SlogP)
    @test get_property(mol, :SlogP) == "1.1878"
end
    
    @testset "Writing" begin
        sys = load_sdfile("data/sdfile_test_1.sdf")
        mols = molecules(sys)

        (single_name, single_file) = mktemp(;cleanup = true)

        write_sdfile(single_name, mols[1])
        m_sd = molecules(load_sdfile(single_name))[1]

    # since the molecules have different systems, we cannot simply compare them directly
    # also, m_sd contains the atom_idx property due to the GraphMol conversion
    m_sd.properties = filter(((k,v),) -> k != :atom_idx, m_sd.properties)
    @test _compare_without_system(m_sd, mols[1])

        (set_name, set_file) = mktemp(;cleanup = true)

        write_sdfile(set_name, sys)
        ms_sd = molecules(load_sdfile(set_name))

    ms_sd = map(m -> begin m.properties = filter(((k,v),) -> k != :atom_idx, m.properties); m end, ms_sd)

        @test all([_compare_without_system(ms_sd[i], mols[i]) for i in 1:length(ms_sd)])
    end
end