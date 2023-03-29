@testitem "Simple Protein" begin
    using DataFrames

    mol = Protein()

    @test mol isa Protein{Float32}
    @test mol isa AbstractMolecule{Float32}
    @test mol.name == ""
    @test atoms_df(mol) isa AbstractDataFrame
    @test size(atoms_df(mol)) == (0, 14)
    @test length(atoms(mol)) == 0
    @test bonds_df(mol) isa AbstractDataFrame
    @test size(bonds_df(mol)) == (0, 5)
    @test length(bonds(mol)) == 0
    @test residues_df(mol) isa AbstractDataFrame
    @test size(residues_df(mol), 1) == 0
    @test length(residues(mol)) == 0
    @test natoms(mol) == 0
    @test nbonds(mol) == 0
    @test length(mol.properties) == 0
    mol.name = "my_fancy_protein"
    @test mol.name == "my_fancy_protein"
end

@testitem "Filled Protein" begin
    using DataFrames
    
    mol = Protein("my_fancy_molecule")

    # add atoms
    r_tmp = Vector3{Float32}(1.0, 2.0, 4.0)
    adf = DataFrame(
        idx = fill(0, 6),
        number = [i for i in 1:6],
        element = fill(Elements.H, 6),
        name = fill("H", 6),
        atomtype = fill("na", 6),
        r = [i .* Vector3{Float32}(1.0, 2.0, 4.0) for i in 1:6],
        v = fill(Vector3{Float32}(0.0, 0.0, 0.0), 6), 
        F = fill(Vector3{Float32}(0.0, 0.0, 0.0), 6),
        formal_charge = 1,
        charge = 2.0f32,
        radius = 1.02f32,
        has_velocity = fill(false, 6),
        has_force = fill(false, 6),
        properties = Properties()
    )

    ca = Chain(mol, "A")
    r1 = Residue(ca, 1, AminoAcid('A'))
    cb = Chain(mol, "B")
    r2 = Residue(cb, 2, AminoAcid('D'))
    for (i, atom) in enumerate(eachrow(adf))
        push!(i < 4 ? r1 : r2, copy(atom)::AtomTuple{Float32})
    end

    # add bonds
    aidx = Dict(a.number => a.idx for a in eachrow(atoms_df(mol)))
    for i in 1:4
        push!(mol, (
            idx = 0,
            a1 = aidx[i],
            a2 = aidx[5],
            order = BondOrder.Single,
            properties = Properties()
        )::BondTuple)
    end

    @test mol.name == "my_fancy_molecule"
    @test nbonds(mol) == 4
    @test natoms(mol) == 6
    @test natoms(r1) == 3
    @test natoms(r2) == 3

    # test properties
    props = mol.properties
    props["molecule computed"] = false
    props["molecule resolution"] = 2.5
    @test length(props) == 2

    props["resolution unit"] = "Angstroem"
    @test length(props) == 3
    @test !props["molecule computed"]
    @test haskey(props, "resolution unit")
end
