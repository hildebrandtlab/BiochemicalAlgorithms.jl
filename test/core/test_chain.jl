@testitem "Chain" begin
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

            Chain(mol_ds; name = "something", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        chain2 = Chain(mol2; name = "something", properties = Properties(:a => 1), flags = Flags([:A, :B]))

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
        
        @test chain.molecule_idx isa Int
        @test chain.molecule_idx == mol.idx

        @test chain2.name == "something"
        @test chain2.properties == Properties(:a => 1)
        @test chain2.flags == Flags([:A, :B])
        @test chain2.molecule_idx == mol2.idx

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
        chain3.molecule_idx = 999
        @test chain3.molecule_idx == 999

        # chain_by_idx
        @test_throws KeyError chain_by_idx(sys, -1)
        @test chain_by_idx(sys, chain.idx) isa Chain{T}
        @test chain_by_idx(sys, chain.idx) == chain
        
        # chains
        cv = chains(sys)
        @test cv isa ChainTable{T}
        @test length(cv) == 2
        @test length(chains(sys)) == 2
        @test length(chains(sys, molecule_idx = -1)) == 0
        @test length(chains(sys, molecule_idx = mol.idx)) == 1
        @test length(chains(sys, molecule_idx = mol2.idx)) == 1
        @test length(chains(sys, molecule_idx = nothing)) == 2

        # nchains + push!
        @test nchains(sys) isa Int
        @test nchains(sys) == 2
        @test nchains(sys, molecule_idx = -1) == 0
        @test nchains(sys, molecule_idx = mol.idx) == 1
        @test nchains(sys, molecule_idx = mol2.idx) == 1
        @test nchains(sys, molecule_idx = nothing) == 2

        @test Chain(mol).molecule_idx == mol.idx
        @test nchains(sys) == 3
        @test nchains(sys, molecule_idx = -1) == 0
        @test nchains(sys, molecule_idx = mol.idx) == 2
        @test nchains(sys, molecule_idx = mol2.idx) == 1
        @test nchains(sys, molecule_idx = nothing) == 3

        # molecule chains
        mol3 = Molecule(sys)
        @test length(chains(mol3)) == 0
        @test chains(mol3) == chains(sys, molecule_idx = mol3.idx)
        @test nchains(mol3) == 0
        @test nchains(mol3) == nchains(sys, molecule_idx = mol3.idx)

        @test Chain(mol3).molecule_idx == mol3.idx
        @test length(chains(mol3)) == 1
        @test chains(mol3) == chains(sys, molecule_idx = mol3.idx)
        @test nchains(mol3) == 1
        @test nchains(mol3) == nchains(sys, molecule_idx = mol3.idx)

        # chain atoms
        @test length(atoms(chain)) == 0
        @test atoms(chain) == atoms(sys, chain_idx = chain.idx)
        @test natoms(chain) == 0
        @test natoms(chain) == natoms(sys, chain_idx = chain.idx)

        frag = Fragment(chain, 1)
        Atom(frag, 1, Elements.H)
        @test length(atoms(chain)) == 1
        @test atoms(chain) == atoms(sys, chain_idx = chain.idx)
        @test natoms(chain) == 1
        @test natoms(chain) == natoms(sys, chain_idx = chain.idx)

        nuc = Nucleotide(chain, 1)
        Atom(nuc, 2, Elements.C)
        @test length(atoms(chain)) == 2
        @test atoms(chain) == atoms(sys, chain_idx = chain.idx)
        @test natoms(chain) == 2
        @test natoms(chain) == natoms(sys, chain_idx = chain.idx)

        res = Residue(chain, 1, AminoAcid('A'))
        Atom(res, 3, Elements.O)
        @test length(atoms(chain)) == 3
        @test atoms(chain) == atoms(sys, chain_idx = chain.idx)
        @test natoms(chain) == 3
        @test natoms(chain) == natoms(sys, chain_idx = chain.idx)

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
        @test length(bonds(chain)) == 0
        @test bonds(chain) == bonds(sys, chain_idx = chain.idx)
        @test nbonds(chain) == 0
        @test nbonds(chain) == nbonds(sys, chain_idx = chain.idx)

        Bond(chain, Atom(frag, 1, Elements.H).idx, Atom(frag, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(chain)) == 1
        @test bonds(chain) == bonds(sys, chain_idx = chain.idx)
        @test nbonds(chain) == 1
        @test nbonds(chain) == nbonds(sys, chain_idx = chain.idx)

        Bond(chain, Atom(nuc, 1, Elements.H).idx, Atom(nuc, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(chain)) == 2
        @test bonds(chain) == bonds(sys, chain_idx = chain.idx)
        @test nbonds(chain) == 2
        @test nbonds(chain) == nbonds(sys, chain_idx = chain.idx)

        Bond(chain, Atom(res, 1, Elements.H).idx, Atom(res, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(chain)) == 3
        @test bonds(chain) == bonds(sys, chain_idx = chain.idx)
        @test nbonds(chain) == 3
        @test nbonds(chain) == nbonds(sys, chain_idx = chain.idx)
    end
end
