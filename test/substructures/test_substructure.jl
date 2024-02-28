@testitem "Substructure" begin
    # read a simple molecule
    mol = molecules(load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json")))[1]

    filter_fn = atom -> atom.element == Elements.C

    # and create a substructure
    filtered_atoms = filter(filter_fn, atoms(mol))
    idxset = Set(filtered_atoms.idx)
    filtered_bonds = filter(bond ->
        bond.a1 ∈ idxset && bond.a2 ∈ idxset,
        bonds(mol)
    )

    s = Substructure(
        "aspirin substructure",
        mol,
        filtered_atoms,
        filtered_bonds
    )

    @test s isa Substructure
    @test s.name == "aspirin substructure"
    @test natoms(s) == 9
    @test nbonds(s) == 8
    @test length(s.properties) == 3

    # generate the same substructure by filtering
    s2 = filter_atoms(
        filter_fn,
        mol;
        name="aspirin substructure",
        adjacent_bonds=false
    )

    @test s2 == s
end
