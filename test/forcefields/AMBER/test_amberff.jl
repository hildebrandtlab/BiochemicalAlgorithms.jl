# This test item verifies the AmberFF implementation against reference values from C++ BALL tests.
# It covers energy and force calculations for rather simple peptide and small systems.
@testitem "AmberFF" begin
    function _rms_F(sys::System{T}) where T
        √(sum(squared_norm.(atoms(sys).F)) / (3 * natoms(sys)))
    end

    @test_logs BiochemicalAlgorithms.CubicSwitchingFunction{Float64}(0.0, 0.0)
    @test_logs BiochemicalAlgorithms.CubicSwitchingFunction{Float64}(15.0, 15.0)
    @test_logs (:warn, r"cuton value should be smaller than cutoff") BiochemicalAlgorithms.CubicSwitchingFunction{Float64}(15.0, 16.0)

    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)

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
        @test _rms_F(p) ≈ 1703.33f0

        # count_warnings: AlaAla produces 2 torsion warnings
        @test count_warnings(a_ff) == 2
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
            label_terminal_fragments!(sys, fdb)

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
            label_terminal_fragments!(sys, fdb)

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
            label_terminal_fragments!(sys, fdb)

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
            label_terminal_fragments!(sys, fdb)

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
            label_terminal_fragments!(sys, fdb)

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
            @test _rms_F(sys) ≈ T(35.035943894791984)
        end

        # Energy test 5 (AlaGlySer) [AMBER96]
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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
            @test _rms_F(sys) ≈ T(40.63998952652902)
        end

        # Energy test 6 (AlaGlySer) [AMBER94]
        let sys = load_hinfile(ball_data_path("../test/data/AlaGlySer2.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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
            @test _rms_F(sys) ≈ T(41.758379137558485)
        end

        # Force test 1 (Torsion) [AMBER94]
        let sys = load_hinfile(ball_data_path("../test/data/AMBER_test_1.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

            ff = AmberFF(sys, ff_amber94; overwrite_typenames = true)

            # remove non-torsion components
            deleteat!(ff.components, 4)
            deleteat!(ff.components, 2)
            deleteat!(ff.components, 1)
            @test only(ff.components) isa TorsionComponent{T}

            compute_forces!(ff)
            @test atom_by_name(sys, "C").F ≈ T[34.077576, 0.042225838, -19.650848] atol = 1e-4
            @test atom_by_name(sys, "O").F ≈ T[-22.687292, 1.4305115f-6, 13.098242] atol = 1e-5
            @test atom_by_name(sys, "N").F ≈ T[5.8381157, 15.915438, 5.462455] atol = 1e-5
            @test atom_by_name(sys, "H").F ≈ T[-17.228397, -15.957666, 1.0901518] atol = 1e-5
        end

        # Force test 2 (ES switching) [AMBER91/CDIEL]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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
                @test force ≈ dE atol = 10
            end
        end

        # Force test 3 (ES switching) [AMBER91/RDIEL]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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
                @test force ≈ dE atol = 3
            end
        end

        # Force test 4 (VdW switching) [AMBER91]
        let sys = load_hinfile(ball_data_path("../test/data/AmberFF_test_3.hin"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

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

# This testitem verifies the AmberFF implementation against reference values from C++ BALL for the BPTI molecule.
# Reference values from C++ BALL (AmberFF, amber96) for test/data/AmberFF_bench.pdb (BPTI, 892 atoms, SSBOND 5-55, 14-38, 30-51).
# BALL version: 1.7.4 (BALL commit 54e7895fc).
# Generated by test/data/AmberFF_bench_BALL_reference.C with the same preprocessing as
# infer_topology! (normalize_names, add_hydrogens/reconstruct, build_bonds).
# Per-atom data:
# test/data/AmberFF_bench_BALL_reference.tsv (forces in kJ/(mol*A), i.e. BALL's Newton values scaled by N_A / 1e13).
#
# NOTE: C++ BALL does not compute any improper torsions with amber96: it looks up "RES ATOM -" in [ResidueImproperTorsions],
# whose keys are "RES:ATOM", so no key ever matches. BALL's torsion energy and forces therefore contain only
# proper torsions, and impropers are excluded from the comparison below.
#
# C++ BALL computes in single precision (Real = float), so Float64 is compared
# against a Float32 reference with correspondingly loose tolerances.

@testitem "AmberFF vs. BALL reference (BPTI disulfides)" begin

    BALL_E = Dict(
        "Bond Stretches" => 494.5741581,
        "Angle Bends" => 517.217225,
        "Proper Torsion" => 1171.522446,   # BALL "Amber Torsion" (propers only, see above)
        "Van der Waals" => -505.8365109,
        "Electrostatic" => -7374.313773,
    )
    BALL_E_TOTAL = -5696.836454
    BALL_RMS_GRADIENT = 62.3162655

    BALL_SSBOND_RESIDUES = [5, 14, 30, 38, 51, 55]
    BALL_SS_BONDS = Dict((5, 55) => 2.050685644, (14, 38) => 2.088554621, (30, 51) => 2.023274183)

    # (res_id, atom name) => (full_name, AMBER type, charge, F)
    ref = Dict(
        (parse(Int, r[2]), String(r[6])) =>
            (full_name=String(r[5]), type=String(r[7]), charge=parse(Float64, r[8]),
                F=Vector3{Float64}(parse.(Float64, r[9:11])...))
        # skip '#' comment lines (BALL version info) and the column header
        for r in (split(l, '\t') for l in Iterators.drop(
            Iterators.filter(l -> !startswith(l, '#'),
                readlines(ball_data_path("../test/data/AmberFF_bench_BALL_reference.tsv"))), 1))
    )

    for T in [Float32, Float64]
        # BALL is single precision; Julia Float32 matches it more tightly than Float64
        e_rtol = T == Float32 ? 5e-4 : 1e-5
        f_atol = T == Float32 ? 0.05 : 0.1

        p = load_pdb(ball_data_path("../test/data/AmberFF_bench.pdb"), T)
        infer_topology!(p, FragmentDB{T}())
        @test natoms(p) == length(ref)

        # --- disulfides: SSBOND flags, variant names and S-S bonds as in BALL
        @test sort([f.number for f in fragments(p) if has_flag(f, :HAS_SSBOND)]) == BALL_SSBOND_RESIDUES
        for f in fragments(p)
            f.number in BALL_SSBOND_RESIDUES && @test get_full_name(f) == "CYS-S"
        end

        idx_to_atom = Dict(a.idx => a for a in atoms(p))
        ss_bonds = Dict{Tuple{Int,Int},T}()
        for b in bonds(p)
            a1, a2 = idx_to_atom[b.atom1_idx], idx_to_atom[b.atom2_idx]
            (a1.element == Elements.S && a2.element == Elements.S) || continue
            r1, r2 = minmax(parent_fragment(a1).number, parent_fragment(a2).number)
            ss_bonds[(r1, r2)] = distance(a1, a2)
        end
        @test sort(collect(keys(ss_bonds))) == sort(collect(keys(BALL_SS_BONDS)))
        for (k, d) in BALL_SS_BONDS
            @test ss_bonds[k] ≈ d rtol = 1e-5
        end

        # --- typing: every atom gets BALL's full residue name, AMBER type and charge
        ff = AmberFF(p)
        for a in atoms(p)
            r = ref[(parent_fragment(a).number, a.name)]
            @test get_full_name(parent_fragment(a)) == r.full_name
            @test a.atom_type == r.type
            @test a.charge ≈ r.charge atol = 1e-6
        end

        # --- energies
        compute_energy!(ff)
        for (name, e) in BALL_E
            @test ff.energy[name] ≈ T(e) rtol = e_rtol
        end
        @test ff.energy["Hydrogen Bonds"] ≈ zero(T)
        @test ff.energy["Improper Torsion"] > zero(T)   # Julia computes impropers, BALL does not
        @test sum(ff.energy[k] for k in keys(BALL_E)) ≈ T(BALL_E_TOTAL) rtol = e_rtol

        # --- forces without improper torsions (BALL has none, see above)
        atoms(p).F .= Ref(zero(Vector3{T}))
        for c in ff.components
            if c isa BiochemicalAlgorithms.TorsionComponent
                foreach(compute_forces!, c.proper_torsions)
            else
                compute_forces!(c)
            end
        end
        for a in atoms(p)
            @test Vector3{Float64}(a.F) ≈ ref[(parent_fragment(a).number, a.name)].F atol = f_atol
        end
        rms = √(sum(squared_norm.(atoms(p).F)) / (3 * natoms(p)))
        @test rms ≈ T(BALL_RMS_GRADIENT) rtol = 1e-5
    end
end
