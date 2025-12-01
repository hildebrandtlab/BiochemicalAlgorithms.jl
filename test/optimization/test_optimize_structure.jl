@testitem "Optimize structure" tags = [:skip_ci] begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_structure!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ -374.3136f0
end

@testitem "Optimize hydrogen positions" tags = [:skip_ci] begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    ff = AmberFF(sys)
    @test compute_energy!(ff) ≈ 1425.5979f0
    @test Int(optimize_hydrogen_positions!(ff).retcode) == 1
    @test compute_energy!(ff) ≈ 1421.212f0
end


@testitem "Optimize structure - test gradients" tags = [:skip_ci] begin
    using FiniteDifferences    
    using LinearAlgebra
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)

    ff = AmberFF(sys)
    r0 = collect(Float64, Iterators.flatten(atoms(ff.system).r))

    # with increasing order, the normed diff is getting worse?
    fd_tech = FiniteDifferences.central_fdm(15, 1)

    finite_grad = grad(fd_tech, x -> BiochemicalAlgorithms._compute_energy_loss!(x, ff), r0)[1]

    # reset system coordinates
    atoms(ff.system).r .= eachcol(reshape(r0, 3, :))

    grad_r0 = zeros(length(r0))
    BiochemicalAlgorithms._compute_grad!(grad_r0, r0, ff)
    analytic_grad = grad_r0

    normed_diff = norm(finite_grad .- analytic_grad) / norm(analytic_grad)
  
    
    #TODO: What tolerances are appropriate here? depending on the order of fd_tevh
    @test finite_grad ≈ analytic_grad rtol=1e-5, atol=1e-8
    @test normed_diff < 1e-5
    
end