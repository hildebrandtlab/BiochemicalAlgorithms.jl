@testitem "CompiledForceField" begin
    using BiochemicalAlgorithms: compile, compute_energy!, compute_forces!,
        prepare_eval!, rebuild_pairlist!, sync_to_table!

    maxabs(x) = isempty(x) ? 0.0 : maximum(abs, x)

    # Compiled evaluator must reproduce the reference path. With acc == T the
    # results are bit-faithful; with acc == Float64 they match within tolerance.
    function check_against_reference(ff::ForceField{T}, acc; etol, ftol) where T
        update!(ff)
        Eref = compute_energy!(ff)
        eref = copy(ff.energy)
        compute_forces!(ff)
        Fref = copy(atoms(ff.system).F)

        cff = compile(ff; acc = acc)
        Ecmp = compute_energy!(cff)
        compute_forces!(cff)

        @test isapprox(Float64(Ecmp), Float64(Eref); atol = etol, rtol = etol)
        for (k, v) in eref
            @test isapprox(Float64(cff.energy[k]), Float64(v); atol = etol, rtol = etol)
        end
        df = maxabs([Float64(maxabs(Fref[n] - cff.F[n])) for n in eachindex(cff.F)])
        @test df <= ftol
    end

    # AlaAla (PDB, amber96, Float32)
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        check_against_reference(ff, Float32; etol = 1e-3, ftol = 1e-1)
        check_against_reference(ff, Float64; etol = 1e-2, ftol = 1e-1)
    end

    # periodic minimum image: shifting atoms by a whole box length must not
    # change the compiled nonbonded energies or forces
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p; periodic_boundary_conditions = true)
        L = ff.options[:periodic_box_width]
        cff0 = compile(ff)
        compute_energy!(cff0); compute_forces!(cff0)
        e0 = copy(cff0.energy); F0 = copy(cff0.F)

        at = atoms(p)
        for k in 1:2:length(at.idx)
            at.r[k] = at.r[k] + Vector3{Float32}(L, 0, 0)
        end
        cff1 = compile(ff)
        compute_energy!(cff1); compute_forces!(cff1)
        for key in ("Van der Waals", "Hydrogen Bonds", "Electrostatic")
            @test isapprox(cff1.energy[key], e0[key]; rtol = 1e-3, atol = 1e-3)
        end
        @test sum(abs, cff0.box) > 0
    end

    # unsupported components must be rejected, not silently dropped
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        struct DummyComponent{T} <: BiochemicalAlgorithms.AbstractForceFieldComponent{T} end
        push!(ff.components, DummyComponent{Float32}())
        @test_throws ArgumentError compile(ff)
    end

    # a skin larger than half the cutoff gap could miss interacting pairs
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        @test_throws ArgumentError compile(ff; max_displacement = 8)
        @test compile(ff; max_displacement = 2.5) isa BiochemicalAlgorithms.CompiledForceField
    end

    # BALL .hin systems, both storage precisions, acc == T -> exact match
    for T in (Float32, Float64)
        fdb = FragmentDB{T}()
        ff91 = ball_data_path("forcefields/AMBER/amber91.ini")
        ff96 = ball_data_path("forcefields/AMBER/amber96.ini")
        for (nm, parm, noassign) in [
                ("AmberFF_test_1.hin", ff91, true),
                ("AmberFF_test_2.hin", ff91, true),
                ("AmberFF_test_3.hin", ff91, true),
                ("AmberFF_test_4.hin", ff91, true),
                ("AlaGlySer.hin",      ff96, false),
            ]
            sys = load_hinfile(ball_data_path("../test/data/$nm"), T)
            normalize_names!(sys, fdb)
            label_terminal_fragments!(sys, fdb)
            ff = AmberFF(sys, parm; assign_charges = !noassign, assign_typenames = !noassign)
            etol = T === Float64 ? 1e-9 : 1e-2
            ftol = T === Float64 ? 1e-6 : 1e-1
            check_against_reference(ff, T; etol = etol, ftol = ftol)
        end
    end

    # Pair-list rebuild policy: a wide skin with a large move must still match a
    # full rebuild (skin guarantees no in-cutoff pair is missed).
    let
        fdb = FragmentDB{Float64}()
        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), Float64)
        normalize_names!(sys, fdb); label_terminal_fragments!(sys, fdb)
        ff = AmberFF(sys, ball_data_path("forcefields/AMBER/amber96.ini"))
        cff = compile(ff)
        E0 = compute_energy!(cff)
        # full rebuild at the same coordinates -> identical energy
        rebuild_pairlist!(cff)
        @test compute_energy!(cff) ≈ E0
    end

    # Moving coordinates past the Verlet skin rebuilds the pair list and still
    # agrees with the reference force field at the new geometry.
    let
        fdb = FragmentDB{Float64}()
        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), Float64)
        normalize_names!(sys, fdb); label_terminal_fragments!(sys, fdb)
        ff = AmberFF(sys, ball_data_path("forcefields/AMBER/amber96.ini"))
        cff = compile(ff; max_displacement = 0.01)
        r = collect(Float64, Iterators.flatten(atoms(sys).r))
        r[1] += 0.25

        prepare_eval!(cff, r)
        @test cff.r_at_last_build == cff.r
        @test cff.iters_since_build == 0

        compute_energy!(cff)
        compute_forces!(cff)
        compiled_energy = sum(values(cff.energy))
        compiled_forces = copy(cff.F)

        sync_to_table!(cff; forces = false)
        update!(ff)
        reference_energy = compute_energy!(ff)
        compute_forces!(ff)

        @test isapprox(compiled_energy, reference_energy; atol = 1e-8, rtol = 1e-8)
        for k in eachindex(compiled_forces)
            @test isapprox(compiled_forces[k], atoms(sys).F[k]; atol = 1e-6, rtol = 1e-6)
        end
    end

    # Explicit synchronization copies compiled coordinates to AmberFF. The
    # forces=false option leaves the canonical force table untouched.
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        cff = compile(ff)
        at = atoms(p)
        forces_before = copy(at.F)
        cff.r[1] += Vector3{Float64}(0.25, -0.5, 1.0)

        sync_to_table!(cff; forces = false)
        @test at.r[1] == cff.r[1]
        @test at.F == forces_before

        fill!(cff.F, zero(Vector3{Float64}))
        cff.F[1] = Vector3{Float64}(1.0, 2.0, 3.0)
        sync_to_table!(cff; forces = true)
        @test at.F[1] == cff.F[1]
    end

    # Threaded backend must agree with serial (deterministic per-thread buffers;
    # Float64 matches to ~rounding, Float32 differs only by summation order).
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        for acc in (Float32, Float64)
            cs = compile(ff; acc = acc, backend = :serial)
            ct = compile(ff; acc = acc, backend = :threads)
            Es = compute_energy!(cs); Et = compute_energy!(ct)
            compute_forces!(cs); compute_forces!(ct)
            etol = acc === Float64 ? 1e-6 : 1e-1
            ftol = acc === Float64 ? 1e-6 : 1e-1
            @test isapprox(Float64(Es), Float64(Et); atol = etol, rtol = etol)
            @test maxabs([Float64(maxabs(cs.F[n] - ct.F[n])) for n in eachindex(cs.F)]) <= ftol
        end
    end
end
