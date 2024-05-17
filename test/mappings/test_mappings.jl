@testitem "Atom_Bijection:" begin
    mol = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
    mol2 = deepcopy(mol) 

    ab = TrivialAtomBijection(mol, mol2)
    @test ab isa TrivialAtomBijection
    @test size(ab.atoms_A) == size(ab.atoms_B)
end

@testitem "Rigid_Mapping: RigidTransform" begin
    using Rotations: RotMatrix3
    
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 2, 3, 4, 5, 6, 7, 8, 9)
    r = RigidTransform(m,v)

    @test m isa Matrix3{Float32}
    @test v isa Vector3{Float32}
    @test r isa RigidTransform
    @test r.rotation isa RotMatrix3{Float32}
    @test r.translation isa Vector3{Float32}
end

@testitem "Rigid_Mapping: Function translate" begin
    s = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
    s2 = deepcopy(s)
    v = Vector3{Float32}(1.0,1.0,1.0)
    translate!(s, v)

    @test natoms(s) == 21

    for i in eachindex(atoms(s))
        @test isapprox(atoms(s).r[i], atoms(s2).r[i]+v)
    end 
end

@testitem "Rigid_Mapping: Function rigid_transform" begin
    sys = System{Float32}()
    for i in 1:3
        atom = Atom(sys, 1, Elements.H; r = Vector3{Float32}(i*1.0, i*2.0, 0.0))
        @test natoms(sys) == i
    end
      
    sys2 = deepcopy(sys) 
    @test natoms(sys) == 3
    @test natoms(sys2) == 3   
 
    # performs counter clockwise rotation by 90 degree with translation by v=(0, 0, 0)
    v = Vector3{Float32}(0, 0, 0)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)
    
    rigid_transform!(sys2, r)
    @test isapprox(atoms(sys2).r[1], Vector3{Float32}(1., 0., -2.), atol=1e-8)
    @test isapprox(atoms(sys2).r[2], Vector3{Float32}(2., 0., -4.), atol=1e-8)
    @test isapprox(atoms(sys2).r[3], Vector3{Float32}(3., 0., -6.), atol=1e-8)

    
   # performs counter clockwise rotation by 90 degree with translation by v=(1, 1, 1)
    sys3 = deepcopy(sys) 
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)
    
    rigid_transform!(sys3, r)
    @test isapprox(atoms(sys3).r[1], Vector3{Float32}(2., 1., -1.), atol=1e-8)
    @test isapprox(atoms(sys3).r[2], Vector3{Float32}(3., 1., -3.), atol=1e-8)
    @test isapprox(atoms(sys3).r[3], Vector3{Float32}(4., 1., -5.), atol=1e-8)
end

@testitem "Rigid_Mapping: Function compute_rmsd" begin
    sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
    sys2 = deepcopy(sys)
    sys3 = deepcopy(sys) 

    @test isapprox(compute_rmsd(sys, sys2),0.0)

    ab = TrivialAtomBijection(sys, sys2)
    @test ab isa TrivialAtomBijection
    @test isapprox(compute_rmsd(ab), 0.0)
  
    translate!(sys3, Vector3{Float32}(1.0, 0, 0))
    ab2 = TrivialAtomBijection(sys,sys3)

    @test isapprox(compute_rmsd(ab2), 1.0)
    @test isapprox(compute_rmsd(sys, sys3),1.0)

    atom = Atom(sys3, 22, Elements.H, r = atoms(sys).r[1])
    @test natoms(sys3) == 22
    @test natoms(sys) == 21
    @test_throws DimensionMismatch compute_rmsd(sys,sys3)
end

@testitem "Rigid_Mapping: Function compute_rmsd_minimizer" begin
    r_A = [[0., 0., 0.], [0., 1., 0.], [3., 1., 0.]]
    sys = System{Float32}()

    # add atoms
    for i in 1:3
        atom = Atom(sys, 1, Elements.H, r = Vector3{Float32}(r_A[i]))
        @test natoms(sys) == i
    end

    # first: consider simple translation of coordinates by a vector v=(1,1,2)
    sys2 = deepcopy(sys) 
    v = Vector3{Float32}(1.0, 1.0, 2.0)
    translate!(sys2, v)


    # check translated system
    for i in eachindex(atoms(sys2))
        @test isapprox(atoms(sys2).r[i], atoms(sys).r[i] + v)
    end

    # compute transformation 
    ab = TrivialAtomBijection{Float32}(sys2, sys)
    
    rt0 = compute_rmsd_minimizer(ab)
    rt1 = compute_rmsd_minimizer(ab, RMSDMinimizerCoutsias)
    rt2 = compute_rmsd_minimizer(ab, RMSDMinimizerKabsch) # calls internally Coutsias

    @test isapprox(rt1.translation, Vector3{Float32}(-1.0, -1.0, -2.0))
    @test isapprox(rt2.translation, Vector3{Float32}(-1.0, -1.0, -2.0))

    #translate back
    # move atoms of sys2 back to positions of sys
    rigid_transform!(sys2,rt1)
    @test isapprox(atoms(sys2).r[1], atoms(sys).r[1], atol=1e-7)
    @test isapprox(atoms(sys2).r[2], atoms(sys).r[2], atol=1e-8)
    @test isapprox(atoms(sys2).r[3], atoms(sys).r[3], atol=1e-8)
  
    # second: consider simple rotation of coordinates without translation
  
    # performs counter clockwise rotation by 90 degree with translation by v=(0, 0, 0)
    v = Vector3{Float32}(0, 0, 0)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)
    
    rigid_transform!(sys2, r)
    
    @test isapprox(atoms(sys2).r[1], Vector3{Float32}(0., 0., 0.), atol=1e-7)
    @test isapprox(atoms(sys2).r[2], Vector3{Float32}(0., 0., -1.), atol=1e-8)
    @test isapprox(atoms(sys2).r[3], Vector3{Float32}(3., 0., -1.), atol=1e-8)


    #third: consider rotation and translation
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)
    sys2 = deepcopy(sys)
    rigid_transform!(sys2,r)
    @test isapprox(atoms(sys2).r[1], Vector3{Float32}(1., 1., 1.), atol=1e-7)
    @test isapprox(atoms(sys2).r[2], Vector3{Float32}(1., 1., 0.), atol=1e-8)
    @test isapprox(atoms(sys2).r[3], Vector3{Float32}(4., 1., 0.), atol=1e-8)

end

@testitem "Rigid_Mapping: Function map_rigid!" begin

    r_A = [[0., 2., 0.], [0., 2., 1.], [3., 2., 1.]]

    sys = System{Float32}()
    # add atoms
    for i in 1:3
        atom = Atom(sys, 1, Elements.H, r = Vector3{Float32}(r_A[i]))
        @test natoms(sys) == i
    end

    sys2 = deepcopy(sys)
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)

    rigid_transform!(sys2, r)
 

    # now we have two system that can be mapped onto each other
    @test isapprox(atoms(sys2).r[1], Vector3{Float32}(1., 1., -1.), atol=1e-7)
    @test isapprox(atoms(sys2).r[2], Vector3{Float32}(1., 2., -1.), atol=1e-8)
    @test isapprox(atoms(sys2).r[3], Vector3{Float32}(4., 2., -1.), atol=1e-8)
    #rmsd before mapping
    @test isapprox(compute_rmsd(sys2,sys),2.081666f0)

    map_rigid!(sys2,sys)
    #after mapping
    @test isapprox(compute_rmsd(sys2,sys),0.942809f0)

    @test isapprox(atoms(sys2).r[1], Vector3{Float32}(0., 2., 1.3333333), atol=1e-7)
    @test isapprox(atoms(sys2).r[2], Vector3{Float32}(0., 2., 0.3333333), atol=1e-8)
    @test isapprox(atoms(sys2).r[3], Vector3{Float32}(3., 2., 0.3333333), atol=1e-8)


    # test mapping only for heavy atoms

    sys3 = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))

    # create a translated and rotated version
    sys4 = deepcopy(sys3)
    v = Vector3{Float32}(2, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    r = RigidTransform(m,v)

    rigid_transform!(sys4, r)

    #rmsd before mapping
    @test isapprox(compute_rmsd(sys4,sys3), 2.4494898f0)

 
    map_rigid!(sys4, sys3; heavy_atoms_only=true)

    #after mapping
    @test isapprox(compute_rmsd(sys4,sys3), 9.472233f-8)

end
