@testitem "TrivialAtomBijection" begin
    for T in [Float32, Float64]
        mol = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)
        mol2 = deepcopy(mol) 

        ab = TrivialAtomBijection(mol, mol2)
        @test ab isa TrivialAtomBijection{T}
        A, B = atoms(ab)
        @test A isa AtomTable{T}
        @test B isa AtomTable{T}
        @test size(A) == size(B)
    end
end

@testitem "RigidTransform" begin
    for T in [Float32, Float64]
        v = Vector3{T}(1, 1, 1)
        m = Matrix3{T}(1, 2, 3, 4, 5, 6, 7, 8, 9)
        r = RigidTransform(m,v)

        @test r isa RigidTransform{T}
        @test r.rotation == m
        @test r.translation == v
    end
end

@testitem "translate!" begin
    for T in [Float32, Float64]
        s = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)
        s2 = deepcopy(s)
        v = Vector3{T}(1, 1, 1)

        @test translate!(s, v) === s
        for (r1, r2) in zip(atoms(s).r, atoms(s2).r)
            @test isapprox(r1, r2 + v)
        end

        at = atoms(s)
        @test translate!(at, v) === at
        for (r1, r2) in zip(at.r, atoms(s2).r)
            @test isapprox(r1, r2 + 2v)
        end
    end
end

@testitem "rigid_transform!" begin
    for T in [Float32, Float64]
        sys = System{T}()
        for i in 1:3
            Atom(sys, 1, Elements.H; r = Vector3{T}(1i, 2i, 0))
        end
        sys2 = deepcopy(sys) 

        # performs counter clockwise rotation by 90 degree with translation by v=(0, 0, 0)
        v = Vector3{T}(0, 0, 0)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m,v)

        rigid_transform!(sys2, r)
        @test isapprox(atoms(sys2).r[1], Vector3{T}(1, 0, -2))
        @test isapprox(atoms(sys2).r[2], Vector3{T}(2, 0, -4))
        @test isapprox(atoms(sys2).r[3], Vector3{T}(3, 0, -6))

        # performs counter clockwise rotation by 90 degree with translation by v=(1, 1, 1)
        sys3 = deepcopy(sys) 
        v = Vector3{T}(1, 1, 1)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m,v)
        
        @test rigid_transform!(sys3, r) === sys3
        @test isapprox(atoms(sys3).r[1], Vector3{T}(2, 1, -1))
        @test isapprox(atoms(sys3).r[2], Vector3{T}(3, 1, -3))
        @test isapprox(atoms(sys3).r[3], Vector3{T}(4, 1, -5))

        at = atoms(deepcopy(sys))
        rigid_transform!(at, r)
        @test isapprox(at.r[1], Vector3{T}(2, 1, -1))
        @test isapprox(at.r[2], Vector3{T}(3, 1, -3))
        @test isapprox(at.r[3], Vector3{T}(4, 1, -5))
    end
end

@testitem "compute_rmsd" begin
    for T in [Float32, Float64]
        sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)
        sys2 = deepcopy(sys)
        sys3 = deepcopy(sys) 

        @test isapprox(compute_rmsd(sys, sys2), 0)
        @test isapprox(compute_rmsd(atoms(sys), atoms(sys2)), 0)
        @test isapprox(compute_rmsd(TrivialAtomBijection(sys, sys2)),  0)
    
        translate!(sys3, Vector3{T}(1, 0, 0))
        @test isapprox(compute_rmsd(sys, sys3), 1)
        @test isapprox(compute_rmsd(atoms(sys), atoms(sys3)), 1)
        @test isapprox(compute_rmsd(TrivialAtomBijection(sys, sys3)), 1)

        Atom(sys3, 22, Elements.H)
        @test_throws DimensionMismatch compute_rmsd(sys,sys3)
    end
end

@testitem "compute_rmsd_minimizer" begin
    using LinearAlgebra: I, norm
    using Statistics: mean

    for T in [Float32, Float64]
        r_A = [[0., 0., 0.], [0., 1., 0.], [3., 1., 0.]]
        r_A .-= Ref(mean(r_A)) # center at origin
        sys = System{T}()

        # add atoms
        for i in 1:3
            atom = Atom(sys, i, Elements.H, r = Vector3{T}(r_A[i]))
        end

        # first: consider simple translation of coordinates by a vector v=(1,1,2)
        sys2 = deepcopy(sys) 
        v = Vector3{T}(1.0, 1.0, 2.0)
        translate!(sys2, v)

        # compute transformation 
        ab = TrivialAtomBijection(sys2, sys)

        rt0 = compute_rmsd_minimizer(ab)
        rt1 = compute_rmsd_minimizer(ab; minimizer = RMSDMinimizerCoutsias)
        rt2 = compute_rmsd_minimizer(ab; minimizer = RMSDMinimizerKabsch) # calls internally Coutsias

        @test rt0 isa RigidTransform{T}
        @test rt1 isa RigidTransform{T}
        @test rt2 isa RigidTransform{T}

        @test isapprox(rt0.rotation, I)
        @test isapprox(rt1.rotation, I)
        @test isapprox(rt2.rotation, I)

        @test isapprox(rt0.translation, -v)
        @test isapprox(rt1.translation, -v)
        @test isapprox(rt2.translation, -v)

        # second: consider simple rotation of coordinates without translation
        # performs counter clockwise rotation by 90 degree with translation by v=(0, 0, 0)
        v = Vector3{T}(0, 0, 0)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m, v)
        
        sys2 = deepcopy(sys)
        rigid_transform!(sys2, r)
        rt0 = compute_rmsd_minimizer(sys2, sys)
        rt1 = compute_rmsd_minimizer(sys2, sys; minimizer = RMSDMinimizerCoutsias)
        rt2 = compute_rmsd_minimizer(sys2, sys; minimizer = RMSDMinimizerKabsch)

        @test isapprox(rt0.rotation, transpose(m))
        @test isapprox(rt1.rotation, transpose(m))
        @test isapprox(rt2.rotation, transpose(m))

        @test abs(norm(rt0.translation)) < 1e-10
        @test abs(norm(rt1.translation)) < 1e-10
        @test abs(norm(rt2.translation)) < 1e-10

        #third: consider rotation and translation
        v = Vector3{T}(1, 1, 1)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m, v)

        sys2 = deepcopy(sys)
        rigid_transform!(sys2, r)
        rt0 = compute_rmsd_minimizer(sys2, sys)
        rt1 = compute_rmsd_minimizer(sys2, sys; minimizer = RMSDMinimizerCoutsias)
        rt2 = compute_rmsd_minimizer(sys2, sys; minimizer = RMSDMinimizerKabsch)

        @test isapprox(rt0.rotation, transpose(m))
        @test isapprox(rt1.rotation, transpose(m))
        @test isapprox(rt2.rotation, transpose(m))

        @test isapprox(rt0.translation, -v)
        @test isapprox(rt1.translation, -v)
        @test isapprox(rt2.translation, -v)
    end
end

@testitem "map_rigid!" begin
    for T in [Float32, Float64]
        r_A = [[0., 2., 0.], [0., 2., 1.], [3., 2., 1.]]
        sys = System{T}()

        # add atoms
        for i in 1:3
            Atom(sys, 1, Elements.H, r = Vector3{T}(r_A[i]))
        end

        sys2 = deepcopy(sys)
        v = Vector3{T}(1, 1, 1)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m,v)

        rigid_transform!(sys2, r)
        @test map_rigid!(sys2, sys) === sys2
        @test compute_rmsd(sys2,sys) < 1e-6

        # test mappings pubchem structure
        sys3 = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"), T)

        # create a translated and rotated version
        sys4 = deepcopy(sys3)
        v = Vector3{T}(2, 1, 1)
        m = Matrix3{T}(1, 0, 0, 0, 0, -1, 0, 1, 0)
        r = RigidTransform(m,v)

        # all atoms
        rigid_transform!(sys4, r)
        @test map_rigid!(sys4, sys3; heavy_atoms_only=false) === sys4
        @test compute_rmsd(sys4,sys3) <= 1e-6

        # heavy atoms only
        at3, at4 = atoms(sys3), atoms(sys4)
        rigid_transform!(at4, r)
        @test map_rigid!(at4, at3; heavy_atoms_only=true) === at4
        @test compute_rmsd(at4, at3) <= 1e-6
    end
end
