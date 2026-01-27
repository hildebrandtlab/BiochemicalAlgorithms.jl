@testitem "AmberFF" begin
    using Unitful
    const prefactor = ustrip(u"kJ/mol/angstrom"/Unitful.Na |> u"N")


    function _rms_F(sys::System{T}) where T
        √(sum(squared_norm.(atoms(sys).F)) / (3 * natoms(sys)))
    end

    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        fdb = FragmentDB()
        normalize_names!(p, fdb)
        reconstruct_fragments!(p, fdb)
        build_bonds!(p, fdb)

        a_ff = AmberFF(p)

        @test compute_energy!(a_ff) ≈ 1425.5979f0

        @test a_ff.energy["Bond Stretches"]   ≈ 1.3630637f0
        @test a_ff.energy["Angle Bends"]      ≈ 5.40766573f0
        @test a_ff.energy["Proper Torsion"]   ≈ 10.7981319f0
        @test a_ff.energy["Improper Torsion"] ≈ 3.99017335f-06
        @test a_ff.energy["Van der Waals"]    ≈ 1493.17578f0
        @test a_ff.energy["Hydrogen Bonds"]   ≈ 0f0
        @test a_ff.energy["Electrostatic"]    ≈ -85.1466827f0

        compute_forces!(a_ff)
        @test _rms_F(p) ≈ 1682.0159f0
    end

    # BALL tests
    # https://github.com/BALL-Project/ball/blob/master/test/AmberFF_test.C
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()
        ff_amber91 = ball_data_path("forcefields/AMBER/amber91.ini")
        ff_amber94 = ball_data_path("forcefields/AMBER/amber94.ini")
        ff_amber96 = ball_data_path("forcefields/AMBER/amber96.ini")

        # Energy test 1 (single stretch) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_1.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(25.552721448)
            @test ff.energy["Bond Stretches"] ≈ T(25.552721448)
            @test ff.energy["Angle Bends"] ≈ zero(T)
            @test ff.energy["Proper Torsion"] ≈ zero(T)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ zero(T)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ zero(T)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(228.72667607023757)
        end

        # Energy test 2 (Bend) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_2.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(7.895900930651099)
            @test ff.energy["Bond Stretches"] ≈ zero(T) atol=1e-9
            @test ff.energy["Angle Bends"] ≈ T(7.895900930651099)
            @test ff.energy["Proper Torsion"] ≈ zero(T)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ zero(T)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ zero(T)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(51.16440255175596)
        end

        # Energy test 3 (VdW) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(21.577213930269455)
            @test ff.energy["Bond Stretches"] ≈ zero(T) atol=1e-9
            @test ff.energy["Angle Bends"] ≈ zero(T) atol=1e-9
            @test ff.energy["Proper Torsion"] ≈ zero(T)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ T(21.577213930269455)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ zero(T)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(67.577074666338)
        end

        # Energy test 4 (Torsion) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_4.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(16.954878332796728)
            @test ff.energy["Bond Stretches"] ≈ T(3.498308129219706e-6)
            @test ff.energy["Angle Bends"] ≈ T(4.618658784515784e-6) atol=1e-6
            @test ff.energy["Proper Torsion"] ≈ T(17.043042650177043)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ T(-0.08817244725585055)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ zero(T)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(14.284901061469098)
        end

        # Energy test 5 (AlaGlySer) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(-314.13298677164846)
            @test ff.energy["Bond Stretches"] ≈ T(3.004505538267881)
            @test ff.energy["Angle Bends"] ≈ T(8.592612636552754)
            @test ff.energy["Proper Torsion"] ≈ T(0.048692818405399324)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ T(23.319381948394568)
            @test ff.energy["Hydrogen Bonds"] ≈ T(-2.289403578006958)
            @test ff.energy["Electrostatic"] ≈ T(-346.8087761352619)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(24.398090094191605)
        end

        # Energy test 5 (AlaGlySer) [AMBER96]
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber96; assign_charges = false, assign_typenames = false)

            @test compute_energy!(ff) ≈ T(128.93387144346946)
            @test ff.energy["Bond Stretches"] ≈ T(3.0037164329529027)
            @test ff.energy["Angle Bends"] ≈ T(8.145414549467489)
            @test ff.energy["Proper Torsion"] ≈ T(8.085844273285822)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ T(52.495307560889245)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ T(57.20358862687402)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(28.500231091492697)
        end

        # Energy test 6 (AlaGlySer) [AMBER94]
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer2.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber94; overwrite_typenames = true)

            @test compute_energy!(ff) ≈ T(-91.22395429471445)
            @test ff.energy["Bond Stretches"] ≈ T(3.754803527450929)
            @test ff.energy["Angle Bends"] ≈ T(9.109164557456147)
            @test ff.energy["Proper Torsion"] ≈ T(14.355668339616349)
            @test ff.energy["Improper Torsion"] ≈ zero(T)
            @test ff.energy["Van der Waals"] ≈ T(45.44361228166956)
            @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
            @test ff.energy["Electrostatic"] ≈ T(-163.88720300090736)

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(29.38050165182444 )
        end

        # Force test 1 (Torsion) [AMBER94]
        let sys = load_hinfile(ball_data_path("../test/data/AMBER_test_1.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber94; overwrite_typenames = true)

            # remove non-torsion components
            deleteat!(ff.components, 4)
            deleteat!(ff.components, 2)
            deleteat!(ff.components, 1)
            @test only(ff.components) isa TorsionComponent{T}

            compute_forces!(ff)
            @test atom_by_name(sys, "C").F ≈ T[34.07756609401984, 0.04222799449598291, -19.650846410165727] atol = 1e-4
            @test atom_by_name(sys, "O").F ≈ T[-22.6872888378592, 5.01346911541134e-7, 13.098240217939106] atol = 1e-5
            @test atom_by_name(sys, "N").F ≈ T[5.8381215669716475, 15.915440114345756, 5.462454461502626] atol = 1e-5
            @test atom_by_name(sys, "H").F ≈ T[-17.228398823132288, -15.957668610188652, 1.0901517307239912] atol = 1e-5
        end

        # Force test 2 (ES switching) [AMBER91/CDIEL]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91;
                assign_charges = false,
                assign_typenames = false,
                vdw_cuton = T(0.05),
                vdw_cutoff = T(0.1),
                electrostatic_cuton = T(2),
                electrostatic_cutoff = T(4),
                distance_dependent_dielectric = false
            )

            a1, a2 = atoms(sys)
            a1.charge = T(1.35)
            a2.charge = T(-1.2)

            for d in 1.5:0.01:4.5
                a2.r = Vector3{T}(d, 0, 0)

                update!(ff)
                compute_forces!(ff)

                force = a2.F[1] 
                dE = compute_energy!(ff)
                a2.r = Vector3{T}(d - 0.0001, 0, 0)
                update!(ff)
                dE = (compute_energy!(ff) - dE) / T(0.0001)
          
                @test force/prefactor ≈ dE atol = 10
            end
        end

        # Force test 3 (ES switching) [AMBER91/RDIEL]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91;
                assign_charges = false,
                assign_typenames = false,
                vdw_cuton = T(0.05),
                vdw_cutoff = T(0.1),
                electrostatic_cuton = T(2),
                electrostatic_cutoff = T(4),
                distance_dependent_dielectric = true
            )

            a1, a2 = atoms(sys)
            a1.charge = T(1.35)
            a2.charge = T(-1.2)

            for d in 1.5:0.01:4.5
                a2.r = Vector3{T}(d, 0, 0)

                update!(ff)
                compute_forces!(ff)

                force = a2.F[1]
                dE = compute_energy!(ff)
                a2.r = Vector3{T}(d - 0.0001, 0, 0)
                update!(ff)
                dE = (compute_energy!(ff) - dE) / T(0.0001)
                @test force/prefactor ≈ dE atol = 3
            end
        end

        # Force test 4 (VdW switching) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)

            ff = AmberFF(sys, ff_amber91;
                assign_charges = false,
                assign_typenames = false,
                vdw_cuton = T(4),
                vdw_cutoff = T(6),
                electrostatic_cuton = T(0.1),
                electrostatic_cutoff = T(0.2)
            )

            a1, a2 = atoms(sys)
            a1.charge = T(1.35)
            a2.charge = T(-1.2)

            for d in 3.5:0.01:6.5
                a2.r = Vector3{T}(d, 0, 0)

                update!(ff)
                compute_forces!(ff)

                force = a2.F[1] 
                dE = compute_energy!(ff)
                a2.r = Vector3{T}(d - 0.0001, 0, 0)
                update!(ff)
                dE = (compute_energy!(ff) - dE) / T(0.0001)
                @test force ≈ dE atol = 0.01
            end
        end 
        
    end
   
end
