@testitem "ChainTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)

        c1 = Chain(mol;
            name = "my chain",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        c2 = Chain(mol)

        ct = chains(sys)

        # AutoHashEquals, copy, and identity
        ct2 = chains(sys)
        @test ct == ct2
        @test isequal(ct, ct2)
        @test hash(ct) == hash(ct2)
        @test ct !== ct2

        ct2 = copy(ct)
        @test ct == ct2
        @test isequal(ct, ct2)
        @test hash(ct) == hash(ct2)
        @test ct !== ct2

        # Tables.jl interface
        @test Tables.istable(typeof(ct))
        @test Tables.columnaccess(typeof(ct))
        @test Tables.schema(ct) isa Tables.Schema
        @test !isnothing(Tables.columns(ct))
        @test !isnothing(Tables.rows(ct))

        # AbstractArray interface
        @test size(ct) == (2, 2)
        @test length(ct) == 2
        @test eltype(ct) == Chain{T}
        @test keys(ct) == [1, 2]

        # getproperty
        ct._sys === sys
        ct._idx == [c1.idx, c2.idx]

        @test ct.idx isa AbstractVector{Int}
        @test ct.idx == [c1.idx, c2.idx]
        @test ct.name isa AbstractVector{String}
        @test ct.name == [c1.name, c2.name]

        @test ct.properties isa AbstractVector{Properties}
        @test ct.properties == [c1.properties, c2.properties]
        @test ct.flags isa AbstractVector{Flags}
        @test ct.flags == [c1.flags, c2.flags]
        @test ct.molecule_idx isa AbstractVector{Int}
        @test ct.molecule_idx == [c1.molecule_idx, c2.molecule_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(ct, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ct, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(ct, :idx) == Tables.getcolumn(ct, 1) == [c1.idx, c2.idx]
        @test Tables.getcolumn(ct, :name) isa AbstractVector{String}
        @test Tables.getcolumn(ct, 2) isa AbstractVector{String}
        @test Tables.getcolumn(ct, :name) == Tables.getcolumn(ct, 2) == [c1.name, c2.name]

        @test Tables.getcolumn(ct, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(ct, :properties) == [c1.properties, c2.properties]
        @test Tables.getcolumn(ct, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(ct, :flags) == [c1.flags, c2.flags]
        @test Tables.getcolumn(ct, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ct, :molecule_idx) == [c1.molecule_idx, c2.molecule_idx]

        # setproperty!
        @test_throws ErrorException ct.idx = [999, 998]
        @test_throws ErrorException ct.name = ["some other", "names"]

        @test_throws ErrorException ct.properties = [Properties(), Properties(:fourth => 997)]
        @test_throws ErrorException ct.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException ct.molecule_idx = [996, 995]

        # getindex
        @test ct[1] === c1
        @test ct[2] === c2
        @test_throws BoundsError ct[0]
        @test_throws BoundsError ct[3]

        # filter
        @test filter(_ -> true, ct) == ct
        @test only(filter(c -> c.idx == c1.idx, ct)) === c1

        # collect
        cv = collect(ct)
        @test cv isa Vector{Chain{T}}
        @test length(cv) == 2

        # atoms
        @test length(atoms(ct)) == 0
        @test natoms(ct) == 0

        a1 = Atom(c1, 1, Elements.H)
        a2 = Atom(c2, 1, Elements.C)
        @test length(atoms(ct)) == 2
        @test natoms(ct) == 2

        # bonds
        @test length(bonds(ct)) == 0
        @test nbonds(ct) == 0

        Bond(sys, a1.idx, a2.idx, BondOrder.Single)
        @test length(bonds(ct)) == 1
        @test nbonds(ct) == 1

        # chains
        @test nchains(ct) == 2

        # fragments
        @test length(fragments(ct)) == 0
        @test nfragments(ct) == 0
        @test length(fragments(ct; variant = FragmentVariant.None)) == 0
        @test nfragments(ct; variant = FragmentVariant.None) == 0
        @test length(nucleotides(ct)) == 0
        @test nnucleotides(ct) == 0
        @test length(residues(ct)) == 0
        @test nresidues(ct) == 0

        Fragment(c1, 1)
        Nucleotide(c1, 1)
        Residue(c2, 1)
        @test length(fragments(ct)) == 3
        @test nfragments(ct) == 3
        @test length(fragments(ct; variant = FragmentVariant.None)) == 1
        @test nfragments(ct; variant = FragmentVariant.None) == 1
        @test length(nucleotides(ct)) == 1
        @test nnucleotides(ct) == 1
        @test length(residues(ct)) == 1
        @test nresidues(ct) == 1
    end
end

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
        @test parent_molecule(chain) === mol

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            parent(chain_ds) === default_system()
            parent_system(chain_ds) === default_system()

            Chain(mol_ds; name = "something", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        chain2 = Chain(mol2; name = "something", properties = Properties(:a => 1), flags = Flags([:A, :B]))

        # getproperty
        @test chain.idx isa Int
        @test chain.name isa String
        @test chain.name == ""

        @test chain.properties isa Properties
        @test chain.properties == Properties()
        @test chain.flags isa Flags
        @test chain.flags == Flags()
        @test chain.molecule_idx isa Int
        @test chain.molecule_idx == mol.idx

        @test chain._sys isa System{T}
        @test chain._idx isa Int

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

        @test push!(mol, chain) === mol
        @test nchains(sys) == 3
        @test nchains(sys, molecule_idx = -1) == 0
        @test nchains(sys, molecule_idx = mol.idx) == 2
        @test nchains(sys, molecule_idx = mol2.idx) == 1
        @test nchains(sys, molecule_idx = nothing) == 3

        lchain = last(chains(mol))
        @test lchain.idx != chain.idx
        @test lchain.name == chain.name
        @test lchain.properties == chain.properties
        @test lchain.flags == chain.flags
        @test lchain.molecule_idx == mol.idx

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

        res = Residue(chain, 1)
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

        # delete!
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 4
        @test nfragments(sys) == 3

        aidx = Set(atoms(chain).idx)
        @test delete!(chain; keep_atoms = true) === nothing
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nfragments(sys) == 0
        @test_throws KeyError chain.idx
        @test all(isnothing, filter(a -> a.idx in aidx, atoms(sys)).chain_idx)

        frag = Fragment(chain2, 1)
        Bond(Atom(frag, 1, Elements.H), Atom(frag, 2, Elements.C), BondOrder.Single)
        @test natoms(sys) == 14
        @test nbonds(sys) == 4
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nfragments(sys) == 1

        @test delete!(chain2; keep_atoms = false) === nothing
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 2
        @test nfragments(sys) == 0
        @test_throws KeyError chain2.idx
    end
end
