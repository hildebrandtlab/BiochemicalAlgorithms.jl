using Statistics: mean
using LinearAlgebra: Hermitian, eigen
@testitem "Atom_Bijection:" begin
    mol = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
    mol2 = deepcopy(mol) 

    ab = TrivialAtomBijection(mol, mol2)
    @test ab isa TrivialAtomBijection
    @test size(ab.atoms_A) == size(ab.atoms_B)

end

@testitem "Rigid_Mapping: RigidTransform" begin
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 2, 3, 4, 5, 6, 7, 8, 9)
    r = RigidTransform(m,v)

    @test r isa RigidTransform
    @test r.rotation isa Matrix3{Float32}
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

    # add atoms
    for i in 1:3
        atom = Atom(sys, 1, Elements.H;
            name = "my fancy atom",
            atom_type = "my atom type",
            r = Vector3{Float32}(i*1.0, i*2.0, 0.0),
            v = ones(Vector3{Float32}),
            F = ones(Vector3{Float32}),
            formal_charge = 4,
            charge = Float32(5.),
            radius = Float32(6.),
            properties = Properties(:first => 'a', :second => "b"),
            flags =  Flags([:third])
        )
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
    
    # check first atom
    @test isapprox(atoms(sys2).r[1][1], 1.0)
    @test isapprox(atoms(sys2).r[1][2], 0.0)
    @test isapprox(atoms(sys2).r[1][3], -2.0)
    
    # check second atom
    @test isapprox(atoms(sys2).r[2][1], 2.0)
    @test isapprox(atoms(sys2).r[2][2], 0.0)
    @test isapprox(atoms(sys2).r[2][3], -4.0)
    
    # check third atom
    @test isapprox(atoms(sys2).r[3][1], 3.0)
    @test isapprox(atoms(sys2).r[3][2], 0.0)
    @test isapprox(atoms(sys2).r[3][3], -6.0)
    
   # performs counter clockwise rotation by 90 degree with translation by v=(1, 1, 1)
    sys3 = deepcopy(sys) 
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)
    
    rigid_transform!(sys3, r)
    
    # check first atom
    @test isapprox(atoms(sys3).r[1][1], 2.0)
    @test isapprox(atoms(sys3).r[1][2], 1.0)
    @test isapprox(atoms(sys3).r[1][3], -1.0)
    
    # check first atom
    @test isapprox(atoms(sys3).r[2][1], 3.0)
    @test isapprox(atoms(sys3).r[2][2], 1.0)
    @test isapprox(atoms(sys3).r[2][3], -3.0)
    
    # check first atom
    @test isapprox(atoms(sys3).r[3][1], 4.0)
    @test isapprox(atoms(sys3).r[3][2], 1.0)
    @test isapprox(atoms(sys3).r[3][3], -5.0) 
        
end

@testitem "Rigid_Mapping: Function compute_rmsd" begin
    sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
    
    sys2 = deepcopy(sys)

    ab = TrivialAtomBijection(sys, sys2)
    @test ab isa TrivialAtomBijection
    @test isapprox(compute_rmsd(ab), 0.0)

    sys3 = deepcopy(sys) 
    translate!(sys3, Vector3{Float32}(1.0, 0, 0))
    ab2 = TrivialAtomBijection(sys,sys3)

    @test isapprox(compute_rmsd(ab2), 1.0)

end

@testitem "Rigid_Mapping: Function compute_rmsd_minimizer" begin
    sys = System{Float32}()

    # add atoms
    for i in 1:3
        atom = Atom(sys, 1, Elements.H;
            name = "my fancy atom",
            atom_type = "my atom type",
            r = Vector3{Float32}(i*1.0, i*2.0, i*3.0),
            v = ones(Vector3{Float32}),
            F = ones(Vector3{Float32}),
            formal_charge = 4,
            charge = Float32(5.),
            radius = Float32(6.),
            properties = Properties(:first => 'a', :second => "b"),
            flags =  Flags([:third])
        )
        @test natoms(sys) == i
    end
      
    sys2 = deepcopy(sys) 
    @test natoms(sys) == 3
    @test natoms(sys2) == 3   

    # first: consider simple translation of coordinates by a vector v=(1,1,0)
    sys2 = deepcopy(sys) 
    v = Vector3{Float32}(1.0, 1.0, 0.0)
    translate!(sys2, v)

    # check first system
    for i in eachindex(atoms(sys))
        @test isapprox(atoms(sys).r[i], Vector3{Float32}(i*1.0, i*2.0, i*3.0))
    end

    # check translated system
    for i in eachindex(atoms(sys2))
        @test isapprox(atoms(sys2).r[i], atoms(sys).r[i] + v)
    end

    # translate back
    ab = TrivialAtomBijection{Float32}(sys2, sys)
    rt = compute_rmsd_minimizer(ab)

    @test rt isa RigidTransform
    @test rt.rotation isa Matrix3{Float32}
    @test rt.translation isa Vector3{Float32}
    @test isapprox(rt.translation, Vector3{Float32}(1.0, 1.0, 0.0))

    #translate back
    # move atoms of sys2 back to positions of sys
    rigid_transform!(sys2,rt)
    @test atoms(sys2).r[1] == Vector3{Float32}(1.0, 2.0, 3.0) 
    @test atoms(sys2).r[2] == Vector3{Float32}(2.0, 4.0, 3.0)
    @test atoms(sys2).r[3] == Vector3{Float32}(3.0, 6.0, 3.0)
    
    # second: consider simple rotation of coordinates without translation
  
 #=  
    r_A = mol.atoms.r
    r_B = mol2.atoms.r
    mean_A = mean(r_A)
    mean_B = mean(r_B)


    R = mapreduce(t -> t[1] * transpose(t[2]), +, zip(r_B .- Ref(mean_B), r_A .- Ref(mean_A)))

    C = Hermitian(transpose(R) * R)
  
    μ, a = eigen(C)

    println("mu ",μ )
    println("a", a)

    ab2 = TrivialAtomBijection{Float64}(mol, mol2)
    rt = compute_rmsd_minimizer(ab2)
    println(rt.translation)
    println(rt.rotation) =#
end
