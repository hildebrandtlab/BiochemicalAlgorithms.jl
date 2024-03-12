@testitem "FragmentTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol)

        f1 = Fragment(chain, 1;
            name = "my fragment",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        f2 = Fragment(chain, 2)

        ft = fragments(sys)

        # Tables.jl interface
        @test Tables.istable(typeof(ft))
        @test Tables.columnaccess(typeof(ft))
        @test Tables.schema(ft) isa Tables.Schema
        @test !isnothing(Tables.columns(ft))
        @test !isnothing(Tables.rows(ft))

        # AbstractArray interface
        @test size(ft) == (2, 3)
        @test length(ft) == 2
        @test eltype(ft) == Fragment{T}
        @test keys(ft) == [1, 2]

        # getproperty
        @test ft._sys === sys
        @test ft._idx == [f1.idx, f2.idx]

        @test ft.idx isa AbstractVector{Int}
        @test ft.idx == [f1.idx, f2.idx]
        @test ft.number isa AbstractVector{Int}
        @test ft.number == [f1.number, f2.number]
        @test ft.name isa AbstractVector{String}
        @test ft.name == [f1.name, f2.name]

        @test ft.properties isa AbstractVector{Properties}
        @test ft.properties == [f1.properties, f2.properties]
        @test ft.flags isa AbstractVector{Flags}
        @test ft.flags == [f1.flags, f2.flags]
        @test ft.molecule_idx isa AbstractVector{Int}
        @test ft.molecule_idx == [f1.molecule_idx, f2.molecule_idx]
        @test ft.chain_idx isa AbstractVector{Int}
        @test ft.chain_idx == [f1.chain_idx, f2.chain_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(ft, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :idx) == Tables.getcolumn(ft, 1) == [f1.idx, f2.idx]
        @test Tables.getcolumn(ft, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :number) == Tables.getcolumn(ft, 2) == [f1.number, f2.number]
        @test Tables.getcolumn(ft, :name) isa AbstractVector{String}
        @test Tables.getcolumn(ft, 3) isa AbstractVector{String}
        @test Tables.getcolumn(ft, :name) == Tables.getcolumn(ft, 3) == [f1.name, f2.name]

        @test Tables.getcolumn(ft, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(ft, :properties) == [f1.properties, f2.properties]
        @test Tables.getcolumn(ft, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(ft, :flags) == [f1.flags, f2.flags]
        @test Tables.getcolumn(ft, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :molecule_idx) == [f1.molecule_idx, f2.molecule_idx]
        @test Tables.getcolumn(ft, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :chain_idx) == [f1.chain_idx, f2.chain_idx]

        # setproperty!
        @test_throws ErrorException ft.idx = [999, 998]
        @test_throws ErrorException ft.number = [997, 996]
        @test_throws ErrorException ft.name = ["some other", "names"]

        @test_throws ErrorException ft.properties = [Properties(), Properties(:fourth => 995)]
        @test_throws ErrorException ft.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException ft.molecule_idx = [994, 993]
        @test_throws ErrorException ft.chain_idx = [992, 991]

        # getindex
        @test ft[1] === f1
        @test ft[2] === f2
        @test_throws BoundsError ft[0]
        @test_throws BoundsError ft[3]

        # filter
        @test filter(_ -> true, ft) == ft
        @test only(filter(f -> f.idx == f1.idx, ft)) === f1

        # collect
        fv = collect(ft)
        @test fv isa Vector{Fragment{T}}
        @test length(fv) == 2
    end
end

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
        @test parent_molecule(frag) === mol
        @test parent_chain(frag) === chain

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
        @test length(getfield(frag, :_row)) == 3

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
        @test frag.molecule_idx isa Int
        @test frag.molecule_idx == mol.idx
        @test frag.chain_idx isa Int
        @test frag.chain_idx == chain.idx

        @test frag._sys isa System{T}
        @test frag._row isa BiochemicalAlgorithms._FragmentTableRow
        
        @test frag2.number == 1
        @test frag2.name == "something"
        @test frag2.properties == Properties(:a => 1)
        @test frag2.flags == Flags([:A, :B])
        @test frag2.molecule_idx == mol2.idx
        @test frag2.chain_idx == chain2.idx

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
        frag3.molecule_idx = 999
        @test frag3.molecule_idx == 999
        frag3.chain_idx = 998
        @test frag3.chain_idx == 998

        # fragment_by_idx
        @test_throws KeyError fragment_by_idx(sys, -1)
        @test fragment_by_idx(sys, frag.idx) isa Fragment{T}
        @test fragment_by_idx(sys, frag.idx) == frag

        # fragments
        fv = fragments(sys)
        @test fv isa FragmentTable{T}
        @test length(fv) == 2
        @test length(fragments(sys)) == 2
        @test length(fragments(sys, molecule_idx = -1)) == 0
        @test length(fragments(sys, molecule_idx = mol.idx)) == 1
        @test length(fragments(sys, molecule_idx = mol2.idx)) == 1
        @test length(fragments(sys, molecule_idx = nothing)) == 2
        @test length(fragments(sys, chain_idx = -1)) == 0
        @test length(fragments(sys, chain_idx = chain.idx)) == 1
        @test length(fragments(sys, chain_idx = chain2.idx)) == 1
        @test length(fragments(sys, chain_idx = nothing)) == 2
        @test length(fragments(sys, molecule_idx = -1, chain_idx = chain.idx)) == 0
        @test length(fragments(sys, molecule_idx = mol.idx, chain_idx = -1)) == 0
        @test length(fragments(sys, molecule_idx = mol.idx, chain_idx = chain.idx)) == 1
        @test length(fragments(sys, molecule_idx = mol.idx, chain_idx = nothing)) == 1
        @test length(fragments(sys, molecule_idx = nothing, chain_idx = chain.idx)) == 1
        @test length(fragments(sys, molecule_idx = nothing, chain_idx = nothing)) == 2

        # nfragments + push!
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 2
        @test nfragments(sys, molecule_idx = -1) == 0
        @test nfragments(sys, molecule_idx = mol.idx) == 1
        @test nfragments(sys, molecule_idx = mol2.idx) == 1
        @test nfragments(sys, molecule_idx = nothing) == 2
        @test nfragments(sys, chain_idx = -1) == 0
        @test nfragments(sys, chain_idx = chain.idx) == 1
        @test nfragments(sys, chain_idx = chain2.idx) == 1
        @test nfragments(sys, chain_idx = nothing) == 2
        @test nfragments(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 1
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = nothing) == 1
        @test nfragments(sys, molecule_idx = nothing, chain_idx = chain.idx) == 1
        @test nfragments(sys, molecule_idx = nothing, chain_idx = nothing) == 2

        @test push!(chain, frag) === chain
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 3
        @test nfragments(sys, molecule_idx = -1) == 0
        @test nfragments(sys, molecule_idx = mol.idx) == 2
        @test nfragments(sys, molecule_idx = mol2.idx) == 1
        @test nfragments(sys, molecule_idx = nothing) == 3
        @test nfragments(sys, chain_idx = -1) == 0
        @test nfragments(sys, chain_idx = chain.idx) == 2
        @test nfragments(sys, chain_idx = chain2.idx) == 1
        @test nfragments(sys, chain_idx = nothing) == 3
        @test nfragments(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 2
        @test nfragments(sys, molecule_idx = mol.idx, chain_idx = nothing) == 2
        @test nfragments(sys, molecule_idx = nothing, chain_idx = chain.idx) == 2
        @test nfragments(sys, molecule_idx = nothing, chain_idx = nothing) == 3

        lfrag = last(fragments(chain))
        @test lfrag.idx != frag.idx
        @test lfrag.number == frag.number
        @test lfrag.name == frag.name
        @test lfrag.properties == frag.properties
        @test lfrag.flags == frag.flags
        @test lfrag.molecule_idx == parent_molecule(chain).idx
        @test lfrag.chain_idx == chain.idx

        # chain/molecule fragments
        mol3 = Molecule(sys)
        @test size(fragments(mol3), 1) == 0
        @test fragments(mol3) == fragments(sys, molecule_idx = mol3.idx)
        @test nfragments(mol3) == 0
        @test nfragments(mol3) == nfragments(sys, molecule_idx = mol3.idx)

        chain3 = Chain(mol3)
        @test size(fragments(chain3), 1) == 0
        @test fragments(chain3) == fragments(sys, chain_idx = chain3.idx)
        @test nfragments(chain3) == 0
        @test nfragments(chain3) == nfragments(sys, chain_idx = chain3.idx)

        Fragment(chain3, 1)
        @test size(fragments(mol3), 1) == 1
        @test fragments(mol3) == fragments(sys, molecule_idx = mol3.idx)
        @test nfragments(mol3) == 1
        @test nfragments(mol3) == nfragments(sys, molecule_idx = mol3.idx)

        @test size(fragments(chain3), 1) == 1
        @test fragments(chain3) == fragments(sys, chain_idx = chain3.idx)
        @test nfragments(chain3) == 1
        @test nfragments(chain3) == nfragments(sys, chain_idx = chain3.idx)

        # fragment atoms
        @test length(atoms(frag)) == 0
        @test atoms(frag) == atoms(sys, fragment_idx = frag.idx)
        @test natoms(frag) == 0
        @test natoms(frag) == natoms(sys, fragment_idx = frag.idx)

        @test Atom(frag, 1, Elements.H).fragment_idx == frag.idx
        @test length(atoms(frag)) == 1
        @test atoms(frag) == atoms(sys, fragment_idx = frag.idx)
        @test natoms(frag) == 1
        @test natoms(frag) == natoms(sys, fragment_idx = frag.idx)

        for atom in atoms(frag)
            @test parent_fragment(atom) === frag
        end

        # fragment bonds
        @test length(bonds(frag)) == 0
        @test bonds(frag) == bonds(sys, fragment_idx = frag.idx)
        @test nbonds(frag) == 0
        @test nbonds(frag) == nbonds(sys, fragment_idx = frag.idx)

        Bond(frag, Atom(frag, 1, Elements.H).idx, Atom(frag, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(frag)) == 1
        @test bonds(frag) == bonds(sys, fragment_idx = frag.idx)
        @test nbonds(frag) == 1
        @test nbonds(frag) == nbonds(sys, fragment_idx = frag.idx)
    end
end
