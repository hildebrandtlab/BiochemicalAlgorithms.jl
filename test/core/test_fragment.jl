@testitem "Fragment" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        frag = Fragment(chain, 1)
        @test frag isa Fragment{T}
        @test parent(frag) === sys
        @test parent_system(frag) === sys

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            frag_ds = Fragment(chain_ds, 1)
            parent(frag_ds) === default_system()
            parent_system(frag_ds) === default_system()

            Fragment(chain_ds, 1;
                name = "something",
                properties = Properties(:a => "b"),
                flags = Flags([:A])
            )
        end

        frag2 = Fragment(chain2, 1;
            name = "something",
            properties = Properties(:a => 1),
            flags = Flags([:A, :B])
        )

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(frag, :_row)) == 5

        # getproperty
        @test frag.idx isa Int
        @test frag.number isa Int
        @test frag.number == 1
        @test frag.name isa String
        @test frag.name == ""
        @test frag.properties isa Properties
        @test frag.properties == Properties()
        @test frag.flags isa Flags
        @test frag.flags == Flags()

        @test frag._sys isa System{T}
        @test frag._row isa BiochemicalAlgorithms._FragmentTableRow
        
        @test frag.molecule_id isa Int
        @test frag.molecule_id == mol.idx
        @test frag.chain_id isa Int
        @test frag.chain_id == chain.idx

        @test frag2.number == 1
        @test frag2.name == "something"
        @test frag2.properties == Properties(:a => 1)
        @test frag2.flags == Flags([:A, :B])
        @test frag2.molecule_id == mol2.idx
        @test frag2.chain_id == chain2.idx

        # setproperty!
        frag.number = 0
        @test frag.number == 0
        frag.name = "something else"
        @test frag.name == "something else"
        frag.properties = Properties(:first => "v1", :second => 99)
        @test length(frag.properties) == 2
        @test frag.properties[:first] == "v1"
        @test frag.properties[:second] == 99
        frag.flags = Flags([:C])
        @test length(frag.flags) == 1
        @test :C in frag.flags

        frag3 = Fragment(Chain(Molecule(System{T}())), 1)
        frag3.molecule_id = 999
        @test frag3.molecule_id == 999
        frag3.chain_id = 998
        @test frag3.chain_id == 998

        # fragment_by_idx
        @test_throws KeyError fragment_by_idx(sys, -1)
        @test fragment_by_idx(sys, frag.idx) isa Fragment{T}
        @test fragment_by_idx(sys, frag.idx) == frag

        # fragments
        fv = fragments(sys)
        @test fv isa FragmentTable{T}
        @test length(fv) == 2
        @test length(fragments(sys)) == 2
        @test length(fragments(sys, molecule_id = -1)) == 0
        @test length(fragments(sys, molecule_id = mol.idx)) == 1
        @test length(fragments(sys, molecule_id = mol2.idx)) == 1
        @test length(fragments(sys, molecule_id = nothing)) == 2
        @test length(fragments(sys, chain_id = -1)) == 0
        @test length(fragments(sys, chain_id = chain.idx)) == 1
        @test length(fragments(sys, chain_id = chain2.idx)) == 1
        @test length(fragments(sys, chain_id = nothing)) == 2
        @test length(fragments(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(fragments(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(fragments(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # nfragments + push!
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 2
        @test nfragments(sys, molecule_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx) == 1
        @test nfragments(sys, molecule_id = mol2.idx) == 1
        @test nfragments(sys, molecule_id = nothing) == 2
        @test nfragments(sys, chain_id = -1) == 0
        @test nfragments(sys, chain_id = chain.idx) == 1
        @test nfragments(sys, chain_id = chain2.idx) == 1
        @test nfragments(sys, chain_id = nothing) == 2
        @test nfragments(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = chain.idx) == 1
        @test nfragments(sys, molecule_id = mol.idx, chain_id = nothing) == 1
        @test nfragments(sys, molecule_id = nothing, chain_id = chain.idx) == 1
        @test nfragments(sys, molecule_id = nothing, chain_id = nothing) == 2

        @test Fragment(chain, 1).chain_id == chain.idx
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 3
        @test nfragments(sys, molecule_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx) == 2
        @test nfragments(sys, molecule_id = mol2.idx) == 1
        @test nfragments(sys, molecule_id = nothing) == 3
        @test nfragments(sys, chain_id = -1) == 0
        @test nfragments(sys, chain_id = chain.idx) == 2
        @test nfragments(sys, chain_id = chain2.idx) == 1
        @test nfragments(sys, chain_id = nothing) == 3
        @test nfragments(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = chain.idx) == 2
        @test nfragments(sys, molecule_id = mol.idx, chain_id = nothing) == 2
        @test nfragments(sys, molecule_id = nothing, chain_id = chain.idx) == 2
        @test nfragments(sys, molecule_id = nothing, chain_id = nothing) == 3

        # chain/molecule fragments
        mol3 = Molecule(sys)
        @test size(fragments(mol3), 1) == 0
        @test fragments(mol3) == fragments(sys, molecule_id = mol3.idx)
        @test nfragments(mol3) == 0
        @test nfragments(mol3) == nfragments(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(fragments(chain3), 1) == 0
        @test fragments(chain3) == fragments(sys, chain_id = chain3.idx)
        @test nfragments(chain3) == 0
        @test nfragments(chain3) == nfragments(sys, chain_id = chain3.idx)

        Fragment(chain3, 1)
        @test size(fragments(mol3), 1) == 1
        @test fragments(mol3) == fragments(sys, molecule_id = mol3.idx)
        @test nfragments(mol3) == 1
        @test nfragments(mol3) == nfragments(sys, molecule_id = mol3.idx)

        @test size(fragments(chain3), 1) == 1
        @test fragments(chain3) == fragments(sys, chain_id = chain3.idx)
        @test nfragments(chain3) == 1
        @test nfragments(chain3) == nfragments(sys, chain_id = chain3.idx)

        # fragment atoms
        @test length(atoms(frag)) == 0
        @test atoms(frag) == atoms(sys, fragment_id = frag.idx)
        @test natoms(frag) == 0
        @test natoms(frag) == natoms(sys, fragment_id = frag.idx)

        @test Atom(frag, 1, Elements.H).fragment_id == frag.idx
        @test length(atoms(frag)) == 1
        @test atoms(frag) == atoms(sys, fragment_id = frag.idx)
        @test natoms(frag) == 1
        @test natoms(frag) == natoms(sys, fragment_id = frag.idx)

        # fragment bonds
        @test length(bonds(frag)) == 0
        @test bonds(frag) == bonds(sys, fragment_id = frag.idx)
        @test nbonds(frag) == 0
        @test nbonds(frag) == nbonds(sys, fragment_id = frag.idx)

        Bond(frag, Atom(frag, 1, Elements.H).idx, Atom(frag, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(frag)) == 1
        @test bonds(frag) == bonds(sys, fragment_id = frag.idx)
        @test nbonds(frag) == 1
        @test nbonds(frag) == nbonds(sys, fragment_id = frag.idx)
    end
end
