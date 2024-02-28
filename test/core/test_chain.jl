@testitem "Chain" begin
    using DataFrames

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)

        # constructors + parent
        chain = Chain(mol)
        @test chain isa Chain{T}
        @test parent(chain) === sys
        @test parent_system(chain) === sys

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            parent(chain_ds) === default_system()
            parent_system(chain_ds) === default_system()

            Chain(mol_ds, "something", Properties(:a => "b"), Flags([:A]))
        end

        chain2 = Chain(mol2, "something", Properties(:a => 1), Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(chain, :_row)) == 4

        # getproperty
        @test chain.idx isa Int
        @test chain.name isa String
        @test chain.name == ""
        @test chain.properties isa Properties
        @test chain.properties == Properties()
        @test chain.flags isa Flags
        @test chain.flags == Flags()

        @test chain._sys isa System{T}
        @test chain._row isa BiochemicalAlgorithms._ChainTableRow
        
        @test chain.molecule_id isa Int
        @test chain.molecule_id == mol.idx

        @test chain2.name == "something"
        @test chain2.properties == Properties(:a => 1)
        @test chain2.flags == Flags([:A, :B])
        @test chain2.molecule_id == mol2.idx

        # setproperty!
        chain.name = "something else"
        @test chain.name == "something else"
        chain.properties = Properties(:first => "v1", :second => 99)
        @test length(chain.properties) == 2
        @test chain.properties[:first] == "v1"
        @test chain.properties[:second] == 99
        chain.flags = Flags([:C])
        @test length(chain.flags) == 1
        @test :C in chain.flags

        chain3 = Chain(Molecule(System{T}()))
        chain3.molecule_id = 999
        @test chain3.molecule_id == 999

        # chain_by_idx
        @test_throws KeyError chain_by_idx(sys, -1)
        @test chain_by_idx(sys, chain.idx) isa Chain{T}
        @test chain_by_idx(sys, chain.idx) == chain
        
        # chains_df
        df = chains_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(ChainTuple)))
        @test copy(df[1, :]) isa ChainTuple
        @test size(chains_df(sys), 1) == 2
        @test size(chains_df(sys, molecule_id = -1), 1) == 0
        @test size(chains_df(sys, molecule_id = mol.idx), 1) == 1
        @test size(chains_df(sys, molecule_id = mol2.idx), 1) == 1
        @test size(chains_df(sys, molecule_id = nothing), 1) == 2

        # chains
        cv = chains(sys)
        @test cv isa ChainTable{T}
        @test length(cv) == 2
        @test length(chains(sys)) == 2
        @test length(chains(sys, molecule_id = -1)) == 0
        @test length(chains(sys, molecule_id = mol.idx)) == 1
        @test length(chains(sys, molecule_id = mol2.idx)) == 1
        @test length(chains(sys, molecule_id = nothing)) == 2

        # nchains + push!
        @test nchains(sys) isa Int
        @test nchains(sys) == 2
        @test nchains(sys, molecule_id = -1) == 0
        @test nchains(sys, molecule_id = mol.idx) == 1
        @test nchains(sys, molecule_id = mol2.idx) == 1
        @test nchains(sys, molecule_id = nothing) == 2

        @test push!(mol, ChainTuple()) === mol
        @test nchains(sys) == 3
        @test nchains(sys, molecule_id = -1) == 0
        @test nchains(sys, molecule_id = mol.idx) == 2
        @test nchains(sys, molecule_id = mol2.idx) == 1
        @test nchains(sys, molecule_id = nothing) == 3

        # molecule chains
        mol3 = Molecule(sys)
        @test size(chains_df(mol3), 1) == 0
        @test chains_df(mol3) == chains_df(sys, molecule_id = mol3.idx)
        @test length(chains(mol3)) == 0
        @test chains(mol3) == chains(sys, molecule_id = mol3.idx)
        @test nchains(mol3) == 0
        @test nchains(mol3) == nchains(sys, molecule_id = mol3.idx)

        push!(mol3, ChainTuple())
        @test size(chains_df(mol3), 1) == 1
        @test size(chains_df(mol3)) == size(chains_df(sys, molecule_id = mol3.idx))
        @test length(chains(mol3)) == 1
        @test chains(mol3) == chains(sys, molecule_id = mol3.idx)
        @test nchains(mol3) == 1
        @test nchains(mol3) == nchains(sys, molecule_id = mol3.idx)

        # chain atoms
        @test size(atoms_df(chain), 1) == 0
        @test atoms_df(chain) == atoms_df(sys, chain_id = chain.idx)
        @test length(atoms(chain)) == 0
        @test atoms(chain) == atoms(sys, chain_id = chain.idx)
        @test natoms(chain) == 0
        @test natoms(chain) == natoms(sys, chain_id = chain.idx)

        frag = Fragment(chain, 1)
        push!(frag, AtomTuple{T}(1, Elements.H))
        @test size(atoms_df(chain), 1) == 1
        @test size(atoms_df(chain)) == size(atoms_df(sys, chain_id = chain.idx))
        @test length(atoms(chain)) == 1
        @test atoms(chain) == atoms(sys, chain_id = chain.idx)
        @test natoms(chain) == 1
        @test natoms(chain) == natoms(sys, chain_id = chain.idx)

        nuc = Nucleotide(chain, 1)
        push!(nuc, AtomTuple{T}(2, Elements.C))
        @test size(atoms_df(chain), 1) == 2
        @test size(atoms_df(chain)) == size(atoms_df(sys, chain_id = chain.idx))
        @test length(atoms(chain)) == 2
        @test atoms(chain) == atoms(sys, chain_id = chain.idx)
        @test natoms(chain) == 2
        @test natoms(chain) == natoms(sys, chain_id = chain.idx)

        res = Residue(chain, 1, AminoAcid('A'))
        push!(res, AtomTuple{T}(3, Elements.O))
        @test size(atoms_df(chain), 1) == 3
        @test size(atoms_df(chain)) == size(atoms_df(sys, chain_id = chain.idx))
        @test length(atoms(chain)) == 3
        @test atoms(chain) == atoms(sys, chain_id = chain.idx)
        @test natoms(chain) == 3
        @test natoms(chain) == natoms(sys, chain_id = chain.idx)

        for atom in atoms(chain)
            @test parent_chain(atom) === chain
        end
        @test parent_chain(frag) === chain
        @test parent_chain(nuc) === chain
        @test parent_chain(res) === chain
        @test parent_chain(Atom(frag, 1, Elements.H)) === chain
        @test parent_chain(Atom(nuc, 2, Elements.C)) === chain
        @test parent_chain(Atom(res, 3, Elements.O)) === chain

        # chain bonds
        @test size(bonds_df(chain), 1) == 0
        @test bonds_df(chain) == bonds_df(sys, chain_id = chain.idx)
        @test length(bonds(chain)) == 0
        @test bonds(chain) == bonds(sys, chain_id = chain.idx)
        @test nbonds(chain) == 0
        @test nbonds(chain) == nbonds(sys, chain_id = chain.idx)

        @test push!(chain, BondTuple(
            Atom(frag, 1, Elements.H).idx,
            Atom(frag, 2, Elements.C).idx,
            BondOrder.Single
        )) === chain
        @test size(bonds_df(chain), 1) == 1
        @test size(bonds_df(chain)) == size(bonds_df(sys, chain_id = chain.idx))
        @test length(bonds(chain)) == 1
        @test bonds(chain) == bonds(sys, chain_id = chain.idx)
        @test nbonds(chain) == 1
        @test nbonds(chain) == nbonds(sys, chain_id = chain.idx)

        @test push!(chain, BondTuple(
            Atom(nuc, 1, Elements.H).idx,
            Atom(nuc, 2, Elements.C).idx,
            BondOrder.Single
        )) === chain
        @test size(bonds_df(chain), 1) == 2
        @test size(bonds_df(chain)) == size(bonds_df(sys, chain_id = chain.idx))
        @test length(bonds(chain)) == 2
        @test bonds(chain) == bonds(sys, chain_id = chain.idx)
        @test nbonds(chain) == 2
        @test nbonds(chain) == nbonds(sys, chain_id = chain.idx)

        @test push!(chain, BondTuple(
            Atom(res, 1, Elements.H).idx,
            Atom(res, 2, Elements.C).idx,
            BondOrder.Single
        )) === chain
        @test size(bonds_df(chain), 1) == 3
        @test size(bonds_df(chain)) == size(bonds_df(sys, chain_id = chain.idx))
        @test length(bonds(chain)) == 3
        @test bonds(chain) == bonds(sys, chain_id = chain.idx)
        @test nbonds(chain) == 3
        @test nbonds(chain) == nbonds(sys, chain_id = chain.idx)
    end
end
