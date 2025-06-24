@testitem "Read SDFile" begin
    for T in [Float32, Float64]
        sys = load_sdfile(ball_data_path("../test/data/sdfile_test_1.sdf"), T)
        mols = molecules(sys)

        @test sys isa System{T}
        @test length(mols) == 11

        @test sum(map(natoms, mols)) == 518
        @test sum(map(nbonds, mols)) == 528

        mol = mols[1]
        @test natoms(mol) == 39
        @test nbonds(mol) == 42

        @test has_property(mol, :metadata)
        met = get_property(mol, :metadata)
        @test met["NAME"] == "Abacavir_sulfate"
        @test met["b_1rotN"] == "6"
        @test met["Weight"] == "286.339"
        @test met["TPSA"] == "101.88"
        @test met["a_acc"] == "4"
        @test met["a_don"] == "3"
        @test met["logP(o/w)"] == "0.40906"
        @test met["SlogP"] == "1.1878"
    end
end

@testitem "Write SDFile" begin
    using DataFrames

    @inline function _compare_without_system(m1::AbstractAtomContainer, m2::AbstractAtomContainer)
        m1.name == m2.name &&
            DataFrame(atoms(m1))  == DataFrame(atoms(m2)) &&
            DataFrame(bonds(m1))  == DataFrame(bonds(m2))
    end

    for T in [Float32, Float64]
        sys = load_sdfile(ball_data_path("../test/data/sdfile_test_1.sdf"), T)
        mols = molecules(sys)

        (single_name, single_file) = mktemp(;cleanup = true)

        write_sdfile(single_name, first(mols))
        m_sd = first(molecules(load_sdfile(single_name, T)))

        (set_name, set_file) = mktemp(; cleanup = true)

        write_sdfile(set_name, sys)
        ms_sd = molecules(load_sdfile(set_name, T))

        @test all([_compare_without_system(ms_sd[i], mols[i]) for i in eachindex(ms_sd)])
    end
end
