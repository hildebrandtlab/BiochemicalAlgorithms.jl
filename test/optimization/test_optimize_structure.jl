@testitem "Optimize structure" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ -374.3136f0

    # pre-minimized structure: re-optimizing should barely change energy
    e_before = compute_energy!(ff)
    @test Int(optimize_structure!(ff).retcode) == 1
    e_after = compute_energy!(ff)
    @test isapprox(e_before, e_after; atol=0.01)
end

@testitem "Optimize structure/Float64" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())

    ff = AmberFF(sys)
    @test compute_energy!(ff) isa Float64
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) < -370.0
end

@testitem "Optimize hydrogen positions" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_hydrogen_positions!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ 1421.212f0
end

# ========== MINIBATCH TESTS ==========

@testitem "Optimize structure/minibatch" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)

    ff = AmberFF(sys)
    cff = compile(ff)

    e_before = compute_energy!(ff)
    @test e_before ≈ 1425.5979f0

    sol = optimize_structure_mini!(cff; epochs=5, batchsize=5)
    @test Int(sol.retcode) == 0

    e_after = compute_energy!(cff.ff)
    @test e_after < e_before
end

@testitem "Optimize structure/minibatch Float64" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())

    ff = AmberFF(sys)
    cff = compile(ff)

    sol = optimize_structure_mini!(cff; epochs=3, batchsize=8)
    @test Int(sol.retcode) == 0

    e_final = compute_energy!(cff.ff)
    @test e_final isa Float64
    @test e_final < 1425.0
    @show e_final
end

# ========== BACKEND CONSISTENCY TESTS ==========

@testitem "Optimize structure/minibatch thread consistency" begin
    using Random

    # Setup: create identical systems with serial and threaded backends
    sys1 = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys1)
    ff1 = AmberFF(sys1)
    cff1 = compile(ff1, backend=:serial)

    sys2 = deepcopy(sys1)
    ff2 = AmberFF(sys2)
    cff2 = compile(ff2, backend=:threads)

    # ===== PART 1: Verify force field components match =====
    @test cff1.stretch.i == cff2.stretch.i
    @test cff1.stretch.j == cff2.stretch.j
    @test cff1.bend.i == cff2.bend.i
    @test cff1.proper.i == cff2.proper.i
    @test cff1.improper.i == cff2.improper.i
    @test cff1.es.i == cff2.es.i

    # ===== PART 2: Verify batch creation matches =====
    Random.seed!(1234)
    ds1 = InteractionDataSet(cff1; batchsize=100)
    Random.seed!(1234)
    ds2 = InteractionDataSet(cff2; batchsize=100)

    b1 = ds1.batches[1]
    b2 = ds2.batches[1]

    @test b1.stretch_range == b2.stretch_range
    @test b1.bend_range == b2.bend_range
    @test b1.torsion_range == b2.torsion_range
    @test b1.improper_range == b2.improper_range
    @test b1.es_range == b2.es_range

    # ===== PART 3: Verify _chunk_range partitioning =====
    for i in 1:length(ds1.batches)
        batch1 = ds1.batches[i]
        batch2 = ds2.batches[i]

        for t in 1:Threads.nthreads()
            sr1 = BiochemicalAlgorithms._chunk_range(batch1.stretch_range, Threads.nthreads(), t)
            sr2 = BiochemicalAlgorithms._chunk_range(batch2.stretch_range, Threads.nthreads(), t)
            @test sr1 == sr2

            br1 = BiochemicalAlgorithms._chunk_range(batch1.bend_range, Threads.nthreads(), t)
            br2 = BiochemicalAlgorithms._chunk_range(batch2.bend_range, Threads.nthreads(), t)
            @test br1 == br2

            tr1 = BiochemicalAlgorithms._chunk_range(batch1.torsion_range, Threads.nthreads(), t)
            tr2 = BiochemicalAlgorithms._chunk_range(batch2.torsion_range, Threads.nthreads(), t)
            @test tr1 == tr2

            ir1 = BiochemicalAlgorithms._chunk_range(batch1.improper_range, Threads.nthreads(), t)
            ir2 = BiochemicalAlgorithms._chunk_range(batch2.improper_range, Threads.nthreads(), t)
            @test ir1 == ir2

            hbr1 = BiochemicalAlgorithms._chunk_range(batch1.hb_range, Threads.nthreads(), t)
            hbr2 = BiochemicalAlgorithms._chunk_range(batch2.hb_range, Threads.nthreads(), t)
            @test hbr1 == hbr2
        end
    end

    # Verify all indices are covered exactly once
    batch = ds1.batches[1]
    nt = Threads.nthreads()
    all_threaded_stretch = Int[]
    for t in 1:nt
        sr = BiochemicalAlgorithms._chunk_range(batch.stretch_range, nt, t)
        append!(all_threaded_stretch, collect(sr))
    end
    serial_stretch = collect(batch.stretch_range)
    @test sort(all_threaded_stretch) == sort(serial_stretch)
    @test length(all_threaded_stretch) == length(serial_stretch)

    # ===== PART 4: Verify initial configuration matches =====
    r1 = collect(Float64, Iterators.flatten(atoms(cff1.ff.system).r))
    r2 = collect(Float64, Iterators.flatten(atoms(cff2.ff.system).r))
    @test r1 == r2

    # ===== PART 5: Verify energy computation matches =====
    p1 = BiochemicalAlgorithms.MiniBatchParams(cff1, ds1.batches, 1)
    p2 = BiochemicalAlgorithms.MiniBatchParams(cff2, ds2.batches, 1)

    e1_serial = BiochemicalAlgorithms._compute_energy_loss(r1, p1)
    e2_serial = BiochemicalAlgorithms._compute_energy_loss(r2, p2)
    @test isapprox(e1_serial, e2_serial; rtol=1e-6)

    # ===== PART 6: Verify optimization produces reproducible results =====
    SEED = 12345

    r0_1 = collect(Float64, Iterators.flatten(atoms(cff1.ff.system).r))
    sol1 = optimize_structure_mini!(cff1; epochs=1, batchsize=100, seed=SEED)

    r0_2 = collect(Float64, Iterators.flatten(atoms(cff2.ff.system).r))
    sol2 = optimize_structure_mini!(cff2; epochs=1, batchsize=100, seed=SEED)

    @test sol1.retcode == sol2.retcode
    @test sol1.stats.iterations == sol2.stats.iterations
    @test isapprox(sol1.u, sol2.u; rtol=1e-2)
    @test isapprox(sol1.objective, sol2.objective; rtol=1e-3)
end
