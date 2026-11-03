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

@testitem "Optimize structure/constrained atoms" begin
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
    infer_topology!(sys)
    ff = AmberFF(sys)
    constrained = findall(a -> a.element != Elements.H, atoms(sys))
    r_before = copy(atoms(sys).r)

    old = [1]
    append!(empty!(ff.constrained_atoms), old)
    optimize_hydrogen_positions!(ff; maxiters = 5)
    @test ff.constrained_atoms == old
    @test atoms(sys).r[constrained] == r_before[constrained]
    @test all(iszero, atoms(sys).F[constrained])
    @test any(!iszero, atoms(sys).F)

    # constraints are restored when the optimization throws
    @test_throws ArgumentError optimize_hydrogen_positions!(ff; max_displacement = 100)
    @test ff.constrained_atoms == old
end

@testitem "Optimize structure/callback state" begin
    using BiochemicalAlgorithms: compile, compute_forces!, prepare_eval!

    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())
    ff = AmberFF(sys)
    reference = compile(ff)

    calls = Ref(0)
    mismatch = Ref(0.0)
    optimize_structure!(ff; maxiters = 10, callback = (state, loss) -> begin
        calls[] += 1
        at = atoms(sys)
        mismatch[] = max(mismatch[], maximum(abs, state.u .- collect(Iterators.flatten(at.r))))
        prepare_eval!(reference, state.u)
        compute_forces!(reference)
        mismatch[] = max(mismatch[], maximum(sqrt(squared_norm(reference.F[k] - at.F[k])) for k in eachindex(at.F)))
        false
    end)
    @test calls[] > 0
    @test mismatch[] < 1e-9
end

@testitem "Optimize structure/custom components" begin
    using BiochemicalAlgorithms: AbstractForceFieldComponent, NonBondedComponent, compile

    # a harmonic restraint on the first atom, implemented only through the public
    # component interface
    struct RestraintComponent{T} <: AbstractForceFieldComponent{T}
        ff::ForceField{T}
        energy::Dict{String, T}
    end
    function BiochemicalAlgorithms.compute_energy!(c::RestraintComponent{T}) where T
        e = squared_norm(atoms(c.ff.system).r[1])
        c.energy["Restraint"] = e
        e
    end
    function BiochemicalAlgorithms.compute_forces!(c::RestraintComponent)
        at = atoms(c.ff.system)
        at.F[1] -= 2 * at.r[1]
        nothing
    end

    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys, FragmentDB{Float64}())
    ff = AmberFF(sys)
    push!(ff.components, RestraintComponent{Float64}(ff, Dict{String, Float64}()))
    @test_throws ArgumentError compile(ff)

    e_before = compute_energy!(ff)
    @test_logs (:warn, r"reference path") match_mode=:any optimize_structure!(ff; maxiters = 20)
    @test compute_energy!(ff) < e_before

    # force fields without nonbonded terms stay on the compiled path
    sys2 = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    infer_topology!(sys2, FragmentDB{Float64}())
    ff2 = AmberFF(sys2)
    filter!(c -> !(c isa NonBondedComponent), ff2.components)
    e2 = compute_energy!(ff2)
    optimize_structure!(ff2; maxiters = 20)
    @test compute_energy!(ff2) < e2
end
