@testset "Substructure" begin

    # read a simple molecule
    mol = load_pubchem_json("data/aspirin_pug.json")[1]

    filter_fn = :element => e -> e == Elements.C

    # and create a substructure
    atoms = filter(filter_fn, mol.atoms, view = true)
    bonds = filter([:a1, :a2] => (a1, a2) -> a1 ∈ atoms.number && a2 ∈ atoms.number, mol.bonds, view = true)

    s = Substructure(
        "aspirin substructure",
        mol,
        atoms,
        bonds
    )

    @test s isa Substructure
    @test s.name == "aspirin substructure"
    @test s.atoms isa SubDataFrame
    @test size(s.atoms) == (9,11)
    @test s.bonds isa SubDataFrame
    @test size(s.bonds) == (8,4)
    @test count_atoms(s) == 9
    @test count_bonds(s) == 8
    @test length(s.properties) == 0

    # generate the same substructure by filtering
    s2 = filter_atoms(
        filter_fn,
        mol;
        name="aspirin substructure",
        adjacent_bonds=false
    )

    @test s2 == s
end
