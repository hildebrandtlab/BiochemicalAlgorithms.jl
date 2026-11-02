@testitem "CharmmFF" begin
    function _rms_F(sys::System{T}) where T
        √(sum(squared_norm.(atoms(sys).F)) / (3 * natoms(sys)))
    end

    # maximum relative deviation between the forces and the finite-difference gradient
    function _max_force_error(ff::ForceField{T}, sys::System{T}; h::T = T(1e-6)) where T
        compute_forces!(ff)
        F = copy(atoms(sys).F)

        max_error = zero(T)
        for (i, atom) in enumerate(atoms(sys)), k in 1:3
            r₀ = atom.r
            δ = Vector3{T}(ntuple(j -> j == k ? h : zero(T), 3))

            atom.r = r₀ + δ; update!(ff); E₊ = compute_energy!(ff)
            atom.r = r₀ - δ; update!(ff); E₋ = compute_energy!(ff)
            atom.r = r₀;     update!(ff)

            max_error = max(max_error, abs(-(E₊ - E₋) / 2h - F[i][k]) / max(one(T), abs(F[i][k])))
        end

        max_error
    end

    # BALL tests
    # https://github.com/BALL-Project/ball/blob/master/test/CharmmFF_test1.C
    # https://github.com/BALL-Project/ball/blob/master/test/CharmmFF_test2.C
    for T in [Float32, Float64]
        fdb = FragmentDB{T}()

        function _load_gly()
            sys = load_pdb(ball_data_path("../test/data/CharmmFF_test_1.pdb"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)
            build_bonds!(sys, fdb)
            sys
        end

        # Energy test 1 (GLY) [EEF1]
        let sys = _load_gly()
            ff = CharmmFF(sys)

            @test natoms(sys) == 8
            @test count_warnings(ff) == 0

            @test compute_energy!(ff) ≈ T(157.4572251) atol=5e-3
            @test ff.energy["Bond Stretches"]   ≈ T(30.83411352)  atol=5e-3
            @test ff.energy["Angle Bends"]      ≈ T(307.4544619)  atol=5e-3
            @test ff.energy["Proper Torsion"]   ≈ T(2.11036776)   atol=5e-3
            @test ff.energy["Improper Torsion"] ≈ T(5.48237888)   atol=5e-3
            @test ff.energy["Van der Waals"]    ≈ T(26.64563664)  atol=5e-3
            @test ff.energy["Electrostatic"]    ≈ T(-53.3215644)  atol=5e-3
            @test ff.energy["Solvation"]        ≈ T(-161.7481263) atol=5e-3

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(219.2263702) atol=5e-3
        end

        # Force test 1 (GLY, stretch only) [EEF1]
        let sys = _load_gly()
            ff = CharmmFF(sys)
            filter!(c -> c isa QuadraticStretchComponent, ff.components)

            @test compute_energy!(ff) ≈ T(30.83411352) atol=5e-3

            compute_forces!(ff)
            @test _rms_F(sys) ≈ T(134.2599576) atol=5e-3

            # gradients in kcal/(mol Å) taken from a CHARMM calculation
            # (version 27b1 with param19_eef1.inp and toph19_eef1.inp)
            charmm_gradients = [
                ( 26.44796,  91.42712, -23.82851),
                ( 14.55307,   1.84307, -19.22059),
                (-14.56795, -89.42944,  62.77158),
                (-16.84492, -12.50059,  -5.06838),
                ( -8.71000,   3.45793, -15.72975),
                ( -1.06629, -17.55152,  -0.53075),
                ( 19.60441,  12.80711,  11.19068),
                (-19.41627,   9.94631,  -9.58428)
            ]

            for (atom, ∇E) in zip(atoms(sys), charmm_gradients)
                @test atom.F ≈ -T(4.184) * Vector3{T}(∇E...) atol=5e-2
            end
        end

        function _two_atoms(; kwargs...)
            sys = load_pdb(ball_data_path("../test/data/CharmmFF_test_2.pdb"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)

            ff = CharmmFF(sys; kwargs...)

            a1, a2 = atoms(sys)
            a1.r = zero(Vector3{T})

            ff, a1, a2
        end

        function _check_scan(ff, a1, a2, reference)
            for (d, F_ref, E_ref) in reference
                a2.r = Vector3{T}(d, 0, 0)
                update!(ff)
                compute_forces!(ff)

                @test a1.F[1] ≈ T(F_ref) rtol=1e-4 atol=1e-4
                @test compute_energy!(ff) ≈ T(E_ref) rtol=1e-4 atol=1e-4
            end
        end

        # Force test 9 (2 atoms, ES switching function) [EEF1/CDIEL]
        _check_scan(
            _two_atoms(
                use_eef1 = false,
                vdw_cutoff = 0.1,
                vdw_cuton = 0.05,
                electrostatic_cuton = 2.0,
                electrostatic_cutoff = 4.0,
                distance_dependent_dielectric = false
            )...,
            [(1.2, 1563.02, -1875.625), (2.0, 562.689, -1125.375), (2.2, 605.039, -1008.73),
             (3.0, 703.071, -468.038), (3.8, 261.778, -27.4272), (4.2, 0.0, 0.0)]
        )

        # Force test 10 (2 atoms, vdW switching function) [EEF1]
        _check_scan(
            _two_atoms(
                use_eef1 = false,
                vdw_cutoff = 6.0,
                vdw_cuton = 4.0,
                electrostatic_cuton = 0.05,
                electrostatic_cutoff = 0.1
            )...,
            [(3.2, -8.81132, 0.658678), (3.8, 0.281582, -0.692204), (4.2, 0.589907, -0.497289),
             (5.0, 0.292368, -0.12261), (5.8, 0.0369267, -0.00354749), (6.2, 0.0, 0.0)]
        )

        # Force test 11 (2 atoms, solvation switching function) [EEF1]
        _check_scan(
            _two_atoms(
                use_eef1 = true,
                vdw_cuton = 0.05,
                vdw_cutoff = 0.1,
                electrostatic_cuton = 0.05,
                electrostatic_cutoff = 0.1,
                solvation_cuton = 2.0,
                solvation_cutoff = 4.0
            )...,
            [(1.2, -20.9295, -70.9525), (2.0, -4.68368, -79.0981), (2.2, -4.03846, -79.967),
             (3.0, -2.31205, -82.4715), (3.8, -0.512456, -83.6283), (4.2, 0.0, -83.68)]
        )
    end

    # united atoms: hydrogens bound to extended atom types are removed (unless requested otherwise)
    let fdb = FragmentDB{Float64}()
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
        infer_topology!(sys, fdb)
        @test natoms(sys) == 23

        CharmmFF(sys; remove_extended_hydrogens = false)
        @test natoms(sys) == 23

        ff = CharmmFF(sys)
        @test natoms(sys) == 15
        @test count_warnings(ff) == 0
        @test length(ff.components[3].proper_torsions) == 7
        @test length(ff.components[4].improper_torsions) == 5

        # analytical forces vs. finite differences (impropers, EEF1, distance dependent dielectric)
        @test _max_force_error(ff, sys) < 1e-5
    end

    # all-atom parameters (PARAM22)
    let fdb = FragmentDB{Float64}()
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
        infer_topology!(sys, fdb)

        # PARAM22 uses CHARMM names for the amide hydrogens
        for atom in atoms(sys)
            if atom.name == "H"
                atom.name = "HN"
            end
        end

        ff = CharmmFF(sys, ball_data_path("forcefields/CHARMM/param22.ini"))

        @test natoms(sys) == 23
        @test isempty(ff.unassigned_atoms)
        @test isfinite(compute_energy!(ff))
        @test ff.energy["Solvation"] == 0
        @test length(ff.components[3].proper_torsions) > 0

        @test _max_force_error(ff, sys) < 1e-5
    end
end
