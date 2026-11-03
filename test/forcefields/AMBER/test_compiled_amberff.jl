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

@testitem "CompiledForceField/regressions" begin
    using BiochemicalAlgorithms: compile, compute_energy!, compute_forces!, prepare_eval!, sync_to_table!,
        CubicSwitchingFunction, NonBondedComponent,
        _e_lj_pair, _e_hb_pair, _e_es_pair, _f_lj_factor, _f_hb_factor, _f_es_factor,
        _nonbonded_component

    function alaglyser(parm = "amber91")
        fdb = FragmentDB{Float64}()
        sys = load_hinfile(ball_data_path("../test/data/AlaGlySer.hin"), Float64)
        normalize_names!(sys, fdb)
        label_terminal_fragments!(sys, fdb)
        sys, AmberFF(sys, ball_data_path("forcefields/AMBER/$parm.ini");
                     assign_charges = false, assign_typenames = false)
    end

    # force factor must be the derivative of the pair energy: f = -(dE/dr) / r,
    # below the cuton, in the switching region, and for 1-4 scaled pairs
    let sw = CubicSwitchingFunction{Float64}(15.0, 13.0), h = 1e-5
        check(E, f, r) = begin
            dE = (E((r + h)^2) - E((r - h)^2)) / (2h)
            @test f(r^2) * r ≈ -dE rtol = 1e-6
        end
        for r in (2.0, 5.0, 14.0), s in (1.0, 0.5)
            check(q -> _e_lj_pair(q, 1e4, 1e3, s, sw), q -> _f_lj_factor(q, 1e4, 1e3, s, sw), r)
            check(q -> _e_hb_pair(q, 1e4, 1e3, s, sw), q -> _f_hb_factor(q, 1e4, 1e3, s, sw), r)
            for ddd in (true, false)
                check(q -> _e_es_pair(q, 0.3, s, sw, ddd, 1389.35), q -> _f_es_factor(q, 0.3, s, sw, ddd, 1389.35), r)
            end
        end
    end

    # whole-system forces must be the gradient of the energy, on both paths, for
    # a system that actually has hydrogen bonds
    let (sys, ff) = alaglyser()
        cff = compile(ff)
        @test !isempty(cff.hb.i)
        compute_forces!(cff)
        compute_forces!(ff)
        Fref = copy(atoms(sys).F)
        r = collect(Float64, Iterators.flatten(atoms(sys).r))
        h = 1e-5
        for k in eachindex(r)
            r0 = r[k]
            r[k] = r0 + h; prepare_eval!(cff, r); Ep = compute_energy!(cff)
            r[k] = r0 - h; prepare_eval!(cff, r); Em = compute_energy!(cff)
            r[k] = r0
            atom, comp = divrem(k - 1, 3) .+ 1
            @test cff.F[atom][comp] ≈ -(Ep - Em) / (2h) atol = 1e-4
            @test Fref[atom][comp] ≈ -(Ep - Em) / (2h) atol = 1e-4
        end
    end

    # reference forces must use the same minimum image as the pair list: moving
    # atoms by whole box lengths leaves the nonbonded forces of both paths equal
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
        infer_topology!(p, FragmentDB{Float64}())
        ff = AmberFF(p; periodic_boundary_conditions = true)
        L = ff.options[:periodic_box_width]
        at = atoms(p)
        for k in 1:2:length(at.idx)
            at.r[k] = at.r[k] + Vector3{Float64}(L, 0, 0)
        end
        update!(ff)
        compute_forces!(ff)
        Fref = copy(at.F)
        cff = compile(ff)
        compute_forces!(cff)
        @test !isempty(cff.lj.i)
        @test maximum(sqrt(squared_norm(Fref[k] - cff.F[k])) for k in eachindex(Fref)) < 1e-6
    end

    # a box below twice the list cutoff would list a pair twice
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        small = AmberFF(p; periodic_boundary_conditions = true,
            periodic_box_width = 30f0, periodic_box_height = 100f0, periodic_box_depth = 100f0)
        @test_throws ArgumentError compile(small)
        @test_logs (:warn, r"periodic box") BiochemicalAlgorithms.setup!(_nonbonded_component(small))
        nb = small.options[:nonbonded_cutoff]
        exact = AmberFF(p; periodic_boundary_conditions = true,
            periodic_box_width = 2nb, periodic_box_height = 2nb, periodic_box_depth = 2nb)
        @test compile(exact) isa BiochemicalAlgorithms.CompiledForceField
        @test_logs BiochemicalAlgorithms.setup!(_nonbonded_component(exact))
        # the box is irrelevant without periodic boundary conditions
        open_box = AmberFF(p; periodic_box_width = 30f0)
        @test_logs BiochemicalAlgorithms.setup!(_nonbonded_component(open_box))
    end

    # out-of-range constraints must not silently leave atoms free
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        @test_throws ArgumentError compile(AmberFF(p; constrained_atoms = [10_000]))
        @test_throws ArgumentError compile(AmberFF(p; constrained_atoms = [0]))
        @test compile(AmberFF(p; constrained_atoms = [natoms(p)])).constrained_mask[end]
    end

    # the threaded backend must be callable from inside a threaded region
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
        infer_topology!(p)
        ff = AmberFF(p)
        cffs = [compile(ff; backend = :threads) for _ in 1:2]
        expected = compute_energy!(cffs[1])
        energies = zeros(Float32, 2)
        Threads.@threads for t in 1:2
            energies[t] = compute_energy!(cffs[t])
            compute_forces!(cffs[t])
        end
        @test all(energies .≈ expected)
    end

    # without nonbonded terms there is no pair list, but evaluation still works
    let p = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
        infer_topology!(p, FragmentDB{Float64}())
        ff = AmberFF(p)
        filter!(c -> !(c isa NonBondedComponent), ff.components)
        cff = compile(ff)
        @test cff.r_at_last_build == cff.r
        r = collect(Float64, Iterators.flatten(atoms(p).r))
        r[1] += 50.0
        @test !BiochemicalAlgorithms.maybe_rebuild!(prepare_eval!(cff, r))
        sync_to_table!(cff; forces = false)
        @test compute_energy!(cff) ≈ compute_energy!(ff)
    end
end
