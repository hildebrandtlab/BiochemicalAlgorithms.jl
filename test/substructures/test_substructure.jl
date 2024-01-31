@testitem "Substructure" begin
    using DataFrames
    using Tables, TableOperations

    # read a simple molecule
    mol = molecules(load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json")))[1]

    filter_fn = atom -> atom.element == Elements.C

    # and create a substructure
    filtered_atoms = Tables.materializer(AtomTable{Float32})(TableOperations.filter(filter_fn, BiochemicalAlgorithms._atoms(mol)))
    idxset = Set(filtered_atoms.idx)
    filtered_bonds = filter(
        [:a1, :a2] => 
            (a1, a2) -> a1 ∈ idxset && 
                        a2 ∈ idxset, BiochemicalAlgorithms._bonds(mol), view = true)

    s = Substructure(
        "aspirin substructure",
        mol,
        filtered_atoms,
        filtered_bonds
    )

    @test s isa Substructure
    @test s.name == "aspirin substructure"
    @test atoms_df(s) isa DataFrame
    @test size(atoms_df(s)) == (9,13)
    @test bonds_df(s) isa SubDataFrame
    @test size(bonds_df(s)) == (8,6)
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
