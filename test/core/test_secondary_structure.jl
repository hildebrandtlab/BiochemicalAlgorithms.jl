@testitem "SecondaryStructureTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)

        c = Chain(mol)

        ss1 = SecondaryStructure(
            c,
            1,
            SecondaryStructureElement.Helix;
            name = "H1",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        ss2 = SecondaryStructure(c, 2, SecondaryStructureElement.Coil; name="C1")

        st = secondary_structures(sys)

        # AutoHashEquals, copy, and identity
        st2 = secondary_structures(sys)
        @test st == st2
        @test isequal(st, st2)
        @test hash(st) == hash(st2)
        @test st !== st2

        st2 = copy(st)
        @test st == st2
        @test isequal(st, st2)
        @test hash(st) == hash(st2)
        @test st !== st2

        # Tables.jl interface
        @test Tables.istable(typeof(st))
        @test Tables.columnaccess(typeof(st))
        @test Tables.schema(st) isa Tables.Schema
        @test !isnothing(Tables.columns(st))
        @test !isnothing(Tables.rows(st))

        # AbstractArray interface
        @test size(st) == (2, 4)
        @test length(st) == 2
        @test eltype(st) == SecondaryStructure{T}
        @test keys(st) == [1, 2]

        # getproperty
        st._sys === sys
        st._idx == [ss1.idx, ss2.idx]

        @test st.idx isa AbstractVector{Int}
        @test st.idx == [ss1.idx, ss2.idx]
        @test st.number isa AbstractVector{Int}
        @test st.number == [ss1.number, ss2.number]
        @test st.type isa AbstractVector{SecondaryStructureType}
        @test st.type == [ss1.type, ss2.type]
        @test st.name isa AbstractVector{String}
        @test st.name == [ss1.name, ss2.name]

        @test st.properties isa AbstractVector{Properties}
        @test st.properties == [ss1.properties, ss2.properties]
        @test st.flags isa AbstractVector{Flags}
        @test st.flags == [ss1.flags, ss2.flags]
        @test st.molecule_idx isa AbstractVector{Int}
        @test st.molecule_idx == [ss1.molecule_idx, ss2.molecule_idx]
        @test st.chain_idx isa AbstractVector{Int}
        @test st.chain_idx == [ss1.chain_idx, ss2.chain_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(st, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(st, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(st, :idx) == Tables.getcolumn(st, 1) == [ss1.idx, ss2.idx]
        @test Tables.getcolumn(st, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(st, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(st, :number) == Tables.getcolumn(st, 2) == [ss1.number, ss2.number]
        @test Tables.getcolumn(st, :type) isa AbstractVector{SecondaryStructureType}
        @test Tables.getcolumn(st, 3) isa AbstractVector{SecondaryStructureType}
        @test Tables.getcolumn(st, :type) == Tables.getcolumn(st, 3) == [ss1.type, ss2.type]
        @test Tables.getcolumn(st, :name) isa AbstractVector{String}
        @test Tables.getcolumn(st, 4) isa AbstractVector{String}
        @test Tables.getcolumn(st, :name) == Tables.getcolumn(st, 4) == [ss1.name, ss2.name]

        @test Tables.getcolumn(st, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(st, :properties) == [ss1.properties, ss2.properties]
        @test Tables.getcolumn(st, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(st, :flags) == [ss1.flags, ss2.flags]
        @test Tables.getcolumn(st, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(st, :molecule_idx) == [ss1.molecule_idx, ss2.molecule_idx]
        @test Tables.getcolumn(st, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(st, :chain_idx) == [ss1.chain_idx, ss2.chain_idx]

        # setproperty!
        @test_throws ErrorException st.idx = [999, 998]
        @test_throws ErrorException st.number = [1000, 1001]
        @test_throws ErrorException st.type = [SecondaryStructureElement.Strand, SecondaryStructureElement.Unknown]
        @test_throws ErrorException st.name = ["some other", "names"]

        @test_throws ErrorException st.properties = [Properties(), Properties(:fourth => 997)]
        @test_throws ErrorException st.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException st.molecule_idx = [996, 995]
        @test_throws ErrorException st.chain_idx = [996, 995]

        # getindex
        @test st[1] === ss1
        @test st[2] === ss2
        @test_throws BoundsError st[0]
        @test_throws BoundsError st[3]

        # filter
        @test filter(_ -> true, st) == st
        @test only(filter(ss -> ss.idx == ss1.idx, st)) === ss1

        # collect
        sv = collect(st)
        @test sv isa Vector{SecondaryStructure{T}}
        @test length(sv) == 2

        # secondary_structures
        @test nsecondary_structures(st) == 2

        # fragments
        @test length(fragments(st)) == 0
        @test nfragments(st) == 0
        @test length(fragments(st; variant = FragmentVariant.None)) == 0
        @test nfragments(st; variant = FragmentVariant.None) == 0
        @test length(nucleotides(st)) == 0
        @test nnucleotides(st) == 0
        @test length(residues(st)) == 0
        @test nresidues(st) == 0

        Fragment(ss1, 1)
        Nucleotide(ss1, 1)
        r = Residue(ss2, 1)
        @test length(fragments(st)) == 3
        @test nfragments(st) == 3
        @test length(fragments(st; variant = FragmentVariant.None)) == 1
        @test nfragments(st; variant = FragmentVariant.None) == 1
        @test length(nucleotides(st)) == 1
        @test nnucleotides(st) == 1
        @test length(residues(st)) == 1
        @test nresidues(st) == 1

        # atoms
        @test length(atoms(st)) == 0
        @test natoms(st) == 0

        a1 = Atom(r, 1, Elements.H)
        a2 = Atom(r, 1, Elements.C)
        @test length(atoms(st)) == 2
        @test natoms(st) == 2

        # bonds
        @test length(bonds(st)) == 0
        @test nbonds(st) == 0

        Bond(sys, a1.idx, a2.idx, BondOrder.Single)
        @test length(bonds(st)) == 1
        @test nbonds(st) == 1
    end
end

@testitem "SecondaryStructure" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        ss = SecondaryStructure(chain, 1, SecondaryStructureElement.Coil; name="C1")

        @test ss isa SecondaryStructure{T}
        @test parent(ss) === sys
        @test parent_system(ss) === sys
        @test parent_molecule(ss) === mol
        @test parent_chain(ss) === chain

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            ss_ds = SecondaryStructure(chain_ds, 1, SecondaryStructureElement.Coil; name="C1")
            parent(ss_ds) === default_system()
            parent_system(ss_ds) === default_system()

            SecondaryStructure(chain_ds, 2, SecondaryStructureElement.Helix; name="H1", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        ss2 = SecondaryStructure(chain2, 1, SecondaryStructureElement.Turn; name="T1", properties = Properties(:a => 1), flags = Flags([:A, :B]))

        # getproperty
        @test ss.idx isa Int
        @test ss.number isa Int
        @test ss.number == 1
        @test ss.type isa SecondaryStructureType
        @test ss.type == SecondaryStructureElement.Coil
        @test ss.name isa String
        @test ss.name == "C1"

        @test ss.properties isa Properties
        @test ss.properties == Properties()
        @test ss.flags isa Flags
        @test ss.flags == Flags()
        @test ss.molecule_idx isa Int
        @test ss.molecule_idx == mol.idx
        @test ss.chain_idx isa Int
        @test ss.chain_idx == chain.idx

        @test ss._sys isa System{T}
        @test ss._idx isa Int

        @test ss2.number == 1
        @test ss2.type == SecondaryStructureElement.Turn
        @test ss2.name == "T1"
        @test ss2.properties == Properties(:a => 1)
        @test ss2.flags == Flags([:A, :B])
        @test ss2.molecule_idx == mol2.idx
        @test ss2.chain_idx == chain2.idx

        # setproperty!
        ss.name = "something else"
        @test ss.name == "something else"

        ss.properties = Properties(:first => "v1", :second => 99)
        @test length(ss.properties) == 2
        @test ss.properties[:first] == "v1"
        @test ss.properties[:second] == 99
        ss.flags = Flags([:C])
        @test length(ss.flags) == 1
        @test :C in ss.flags

        ss3 = SecondaryStructure(Chain(Molecule(System{T}())), 1, SecondaryStructureElement.Unknown; name="U1")
        ss3.molecule_idx = 999
        @test ss3.molecule_idx == 999

        # secondary_structure_by_idx
        @test_throws KeyError secondary_structure_by_idx(sys, -1)
        @test secondary_structure_by_idx(sys, ss.idx) isa SecondaryStructure{T}
        @test secondary_structure_by_idx(sys, ss.idx) == ss

        # secondary_structures
        sv = secondary_structures(sys)
        @test sv isa SecondaryStructureTable{T}
        @test length(sv) == 2
        @test length(secondary_structures(sys)) == 2
        @test length(secondary_structures(sys, molecule_idx = -1)) == 0
        @test length(secondary_structures(sys, molecule_idx = mol.idx)) == 1
        @test length(secondary_structures(sys, molecule_idx = mol2.idx)) == 1
        @test length(secondary_structures(sys, molecule_idx = nothing)) == 2

        @test length(secondary_structures(sys, chain_idx = -1)) == 0
        @test length(secondary_structures(sys, chain_idx = chain.idx)) == 1
        @test length(secondary_structures(sys, chain_idx = chain2.idx)) == 1
        @test length(secondary_structures(sys, chain_idx = nothing)) == 2


        # nsecondary_structures + push!
        @test nsecondary_structures(sys) isa Int
        @test nsecondary_structures(sys) == 2
        @test nsecondary_structures(sys, molecule_idx = -1) == 0
        @test nsecondary_structures(sys, molecule_idx = mol.idx) == 1
        @test nsecondary_structures(sys, molecule_idx = mol2.idx) == 1
        @test nsecondary_structures(sys, molecule_idx = nothing) == 2

        @test nsecondary_structures(sys, chain_idx = -1) == 0
        @test nsecondary_structures(sys, chain_idx = chain.idx) == 1
        @test nsecondary_structures(sys, chain_idx = chain2.idx) == 1
        @test nsecondary_structures(sys, chain_idx = nothing) == 2

        @test push!(chain, ss) === chain
        @test nsecondary_structures(sys) == 3
        @test nsecondary_structures(sys, chain_idx = -1) == 0
        @test nsecondary_structures(sys, chain_idx = chain.idx) == 2
        @test nsecondary_structures(sys, chain_idx = chain2.idx) == 1
        @test nsecondary_structures(sys, chain_idx = nothing) == 3

        lss = last(secondary_structures(mol))
        @test lss.idx != ss.idx
        @test lss.number == ss.number
        @test lss.type == ss.type
        @test lss.name == ss.name
        @test lss.properties == ss.properties
        @test lss.flags == ss.flags
        @test lss.molecule_idx == mol.idx
        @test lss.chain_idx == chain.idx

        # molecule secondary structures
        mol3 = Molecule(sys)
        @test length(secondary_structures(mol3)) == 0
        @test secondary_structures(mol3) == secondary_structures(sys, molecule_idx = mol3.idx)
        @test nsecondary_structures(mol3) == 0
        @test nsecondary_structures(mol3) == nsecondary_structures(sys, molecule_idx = mol3.idx)

        @test SecondaryStructure(Chain(mol3), 1, SecondaryStructureElement.Unknown; name="U1").molecule_idx == mol3.idx
        @test length(secondary_structures(mol3)) == 1
        @test secondary_structures(mol3) == secondary_structures(sys, molecule_idx = mol3.idx)
        @test nsecondary_structures(mol3) == 1
        @test nsecondary_structures(mol3) == nsecondary_structures(sys, molecule_idx = mol3.idx)

        # chain atoms
        @test length(atoms(ss)) == 0
        @test natoms(ss) == 0

        frag = Fragment(ss, 1)
        Atom(frag, 1, Elements.H)
        @test length(atoms(ss)) == 1
        @test natoms(ss) == 1

        nuc = Nucleotide(ss, 1)
        Atom(nuc, 2, Elements.C)
        @test length(atoms(ss)) == 2
        @test natoms(ss) == 2

        res = Residue(ss, 1)
        Atom(res, 3, Elements.O)
        @test length(atoms(ss)) == 3
        @test natoms(chain) == 3

        for atom in atoms(ss)
            @test parent_secondary_structure(atom) === ss
        end
        @test parent_secondary_structure(frag) === ss
        @test parent_secondary_structure(nuc) === ss
        @test parent_secondary_structure(res) === ss
        @test parent_secondary_structure(Atom(frag, 1, Elements.H)) === ss
        @test parent_secondary_structure(Atom(nuc, 2, Elements.C)) === ss
        @test parent_secondary_structure(Atom(res, 3, Elements.O)) === ss

        # secondary structure bonds
        @test length(bonds(ss)) == 0
        @test nbonds(ss) == 0

        Bond(chain, Atom(frag, 1, Elements.H).idx, Atom(frag, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(ss)) == 1
        @test nbonds(ss) == 1

        Bond(chain, Atom(nuc, 1, Elements.H).idx, Atom(nuc, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(ss)) == 2
        @test nbonds(ss) == 2

        Bond(chain, Atom(res, 1, Elements.H).idx, Atom(res, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(ss)) == 3
        @test nbonds(ss) == 3

        # delete!
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nsecondary_structures(sys) == 4
        @test nfragments(sys) == 3

        aidx = Set(atoms(ss).idx)
        @test delete!(ss; keep_fragments = true) === nothing
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nsecondary_structures(sys) == 3
        @test nfragments(sys) == 3
        @test_throws KeyError ss.idx
        @test all(isnothing, parent_secondary_structure.(atom_by_idx.(Ref(sys), aidx)))

        frag = Fragment(ss2, 1)
        Bond(Atom(frag, 1, Elements.H), Atom(frag, 2, Elements.C), BondOrder.Single)
        @test natoms(sys) == 14
        @test nbonds(sys) == 4
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nsecondary_structures(sys) == 3
        @test nfragments(sys) == 4

        @test delete!(ss2; keep_fragments = false) === nothing
        @test natoms(sys) == 12
        @test nbonds(sys) == 3
        @test nmolecules(sys) == 3
        @test nchains(sys) == 3
        @test nsecondary_structures(sys) == 2
        @test nfragments(sys) == 3
        @test_throws KeyError ss2.idx
    end
end
