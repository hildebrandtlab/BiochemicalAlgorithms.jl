@testitem "MoleculeTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        m1 = Molecule(sys;
            name = "my molecule",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        p1 = Protein(sys;
            name = "my protein"
        )

        mt = molecules(sys)

        # AutoHashEquals and identity
        mt2 = molecules(sys)
        @test mt == mt2
        @test isequal(mt, mt2)
        @test hash(mt) == hash(mt2)
        @test mt !== mt2

        pt = proteins(sys)
        pt2 = proteins(sys)
        @test pt == pt2
        @test isequal(pt, pt2)
        @test hash(pt) == hash(pt2)
        @test pt !== pt2

        # Tables.jl interface
        @test Tables.istable(typeof(mt))
        @test Tables.columnaccess(typeof(mt))
        @test Tables.schema(mt) isa Tables.Schema
        @test !isnothing(Tables.columns(mt))
        @test !isnothing(Tables.rows(mt))

        # AbstractArray interface
        @test size(mt) == (2, 2)
        @test length(mt) == 2
        @test eltype(mt) == Molecule{T}
        @test keys(mt) == [1, 2]

        # getproperty
        @test mt._sys === sys
        @test mt._idx == [m1.idx, p1.idx]

        @test mt.idx isa AbstractVector{Int}
        @test mt.idx == [m1.idx, p1.idx]
        @test mt.name isa AbstractVector{String}
        @test mt.name == [m1.name, p1.name]

        @test mt.properties isa AbstractVector{Properties}
        @test mt.properties == [m1.properties, p1.properties]
        @test mt.flags isa AbstractVector{Flags}
        @test mt.flags == [m1.flags, p1.flags]
        @test mt.variant isa AbstractVector{MoleculeVariantType}
        @test mt.variant == [m1.variant, p1.variant]

        # Tables.getcolumn
        @test Tables.getcolumn(mt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, :idx) == Tables.getcolumn(mt, 1) == [m1.idx, p1.idx]
        @test Tables.getcolumn(mt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(mt, 2) isa AbstractVector{String}
        @test Tables.getcolumn(mt, :name) == Tables.getcolumn(mt, 2) == [m1.name, p1.name]

        @test Tables.getcolumn(mt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(mt, :properties) == [m1.properties, p1.properties]
        @test Tables.getcolumn(mt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(mt, :flags) == [m1.flags, p1.flags]
        @test Tables.getcolumn(mt, :variant) isa AbstractVector{MoleculeVariantType}
        @test Tables.getcolumn(mt, :variant) == [m1.variant, p1.variant]

        # setproperty!
        @test_throws ErrorException mt.idx = [999, 998]
        @test_throws ErrorException mt.name = ["some other", "name"]

        @test_throws ErrorException mt.properties = [Properties(), Properties(:fourth => 989)]
        @test_throws ErrorException mt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException mt.variant = [MoleculeVariant.Protein, MoleculeVariant.None]

        # getindex
        @test mt[1] === m1
        @test mt[2] === p1
        @test_throws BoundsError mt[0]
        @test_throws BoundsError mt[3]

        # filter
        @test filter(_ -> true, mt) == mt
        @test only(filter(m -> m.idx == m1.idx, mt)) === m1

        # collect
        mv = collect(mt)
        @test mv isa Vector{Molecule{T}}
        @test length(mv) == 2

        # atoms
        @test length(atoms(mt)) == 0
        @test natoms(mt) == 0

        a1 = Atom(m1, 1, Elements.C)
        a2 = Atom(p1, 1, Elements.O)
        @test length(atoms(mt)) == 2
        @test natoms(mt) == 2

        # bonds
        @test length(bonds(mt)) == 0
        @test nbonds(mt) == 0

        Bond(sys, a1.idx, a2.idx, BondOrder.Single)
        @test length(bonds(mt)) == 1
        @test nbonds(mt) == 1

        # molecules
        @test nmolecules(mt) == 2
        @test nmolecules(mt; variant = MoleculeVariant.None) == 1
        @test nproteins(mt) == 1

        # chains
        @test length(chains(mt)) == 0
        @test nchains(mt) == 0

        c1 = Chain(m1)
        c2 = Chain(p1)
        @test length(chains(mt)) == 2
        @test nchains(mt) == 2

        # fragments
        @test length(fragments(mt)) == 0
        @test nfragments(mt) == 0
        @test length(fragments(mt; variant = FragmentVariant.None)) == 0
        @test nfragments(mt; variant = FragmentVariant.None) == 0
        @test length(nucleotides(mt)) == 0
        @test nnucleotides(mt) == 0
        @test length(residues(mt)) == 0
        @test nresidues(mt) == 0

        Fragment(c1, 1)
        Nucleotide(c1, 1)
        Residue(c2, 1)
        @test length(fragments(mt)) == 3
        @test nfragments(mt) == 3
        @test length(fragments(mt; variant = FragmentVariant.None)) == 1
        @test nfragments(mt; variant = FragmentVariant.None) == 1
        @test length(nucleotides(mt)) == 1
        @test nnucleotides(mt) == 1
        @test length(residues(mt)) == 1
        @test nresidues(mt) == 1
    end
end

@testitem "MoleculeTable/None" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        m1 = Molecule(sys;
            name = "my molecule",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        m2 = Molecule(sys)

        # decoy molecule
        Protein(sys)

        mt = molecules(sys; variant = MoleculeVariant.None)

        # Tables.jl interface
        @test Tables.istable(typeof(mt))
        @test Tables.columnaccess(typeof(mt))
        @test Tables.schema(mt) isa Tables.Schema
        @test !isnothing(Tables.columns(mt))
        @test !isnothing(Tables.rows(mt))

        # AbstractArray interface
        @test size(mt) == (2, 2)
        @test length(mt) == 2
        @test eltype(mt) == Molecule{T}
        @test keys(mt) == [1, 2]
        @test all(mol -> mol.variant === MoleculeVariant.None, mt)

        # getproperty
        @test mt._sys === sys
        @test mt._idx == [m1.idx, m2.idx]

        @test mt.idx isa AbstractVector{Int}
        @test mt.idx == [m1.idx, m2.idx]
        @test mt.name isa AbstractVector{String}
        @test mt.name == [m1.name, m2.name]

        @test mt.properties isa AbstractVector{Properties}
        @test mt.properties == [m1.properties, m2.properties]
        @test mt.flags isa AbstractVector{Flags}
        @test mt.flags == [m1.flags, m2.flags]
        @test mt.variant isa AbstractVector{MoleculeVariantType}
        @test mt.variant == [m1.variant, m2.variant]

        # Tables.getcolumn
        @test Tables.getcolumn(mt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(mt, :idx) == Tables.getcolumn(mt, 1) == [m1.idx, m2.idx]
        @test Tables.getcolumn(mt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(mt, 2) isa AbstractVector{String}
        @test Tables.getcolumn(mt, :name) == Tables.getcolumn(mt, 2) == [m1.name, m2.name]

        @test Tables.getcolumn(mt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(mt, :properties) == [m1.properties, m2.properties]
        @test Tables.getcolumn(mt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(mt, :flags) == [m1.flags, m2.flags]
        @test Tables.getcolumn(mt, :variant) isa AbstractVector{MoleculeVariantType}
        @test Tables.getcolumn(mt, :variant) == [m1.variant, m2.variant]

        # setproperty!
        @test_throws ErrorException mt.idx = [999, 998]
        @test_throws ErrorException mt.name = ["some other", "name"]

        @test_throws ErrorException mt.properties = [Properties(), Properties(:fourth => 989)]
        @test_throws ErrorException mt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException mt.variant = [MoleculeVariant.Protein, MoleculeVariant.Protein]

        # getindex
        @test mt[1] === m1
        @test mt[2] === m2
        @test_throws BoundsError mt[0]
        @test_throws BoundsError mt[3]

        # filter
        @test filter(_ -> true, mt) == mt
        @test only(filter(m -> m.idx == m1.idx, mt)) === m1

        # collect
        mv = collect(mt)
        @test mv isa Vector{Molecule{T}}
        @test length(mv) == 2
    end
end

@testitem "MoleculeTable/Protein" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()

        p1 = Protein(sys;
            name = "my protein",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        p2 = Protein(sys)

        # decoy molecule
        Molecule(sys)

        pt = proteins(sys)

        # Tables.jl interface
        @test Tables.istable(typeof(pt))
        @test Tables.columnaccess(typeof(pt))
        @test Tables.schema(pt) isa Tables.Schema
        @test !isnothing(Tables.columns(pt))
        @test !isnothing(Tables.rows(pt))

        # AbstractArray interface
        @test size(pt) == (2, 2)
        @test length(pt) == 2
        @test eltype(pt) == Molecule{T}
        @test keys(pt) == [1, 2]
        @test all(isprotein, pt)

        # getproperty
        @test pt._sys === sys
        @test pt._idx == [p1.idx, p2.idx]

        @test pt.idx isa AbstractVector{Int}
        @test pt.idx == [p1.idx, p2.idx]
        @test pt.name isa AbstractVector{String}
        @test pt.name == [p1.name, p2.name]

        @test pt.properties isa AbstractVector{Properties}
        @test pt.properties == [p1.properties, p2.properties]
        @test pt.flags isa AbstractVector{Flags}
        @test pt.flags == [p1.flags, p2.flags]
        @test pt.variant isa AbstractVector{MoleculeVariantType}
        @test pt.variant == [p1.variant, p2.variant]

        # Tables.getcolumn
        @test Tables.getcolumn(pt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(pt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(pt, :idx) == Tables.getcolumn(pt, 1) == [p1.idx, p2.idx]
        @test Tables.getcolumn(pt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(pt, 2) isa AbstractVector{String}
        @test Tables.getcolumn(pt, :name) == Tables.getcolumn(pt, 2) == [p1.name, p2.name]

        @test Tables.getcolumn(pt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(pt, :properties) == [p1.properties, p2.properties]
        @test Tables.getcolumn(pt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(pt, :flags) == [p1.flags, p2.flags]
        @test Tables.getcolumn(pt, :variant) isa AbstractVector{MoleculeVariantType}
        @test Tables.getcolumn(pt, :variant) == [p1.variant, p2.variant]

        # setproperty!
        @test_throws ErrorException pt.idx = [999, 998]
        @test_throws ErrorException pt.name = ["some other", "name"]

        @test_throws ErrorException pt.properties = [Properties(), Properties(:fourth => 989)]
        @test_throws ErrorException pt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException pt.variant = [MoleculeVariant.None, MoleculeVariant.None]

        # getindex
        @test pt[1] === p1
        @test pt[2] === p2
        @test_throws BoundsError pt[0]
        @test_throws BoundsError pt[3]

        # filter
        @test filter(_ -> true, pt) == pt
        @test only(filter(m -> m.idx == p1.idx, pt)) === p1

        # collect
        mv = collect(pt)
        @test mv isa Vector{Molecule{T}}
        @test length(mv) == 2
    end
end

@testitem "Molecule" begin
    for T in [Float32, Float64]
        sys = System{T}()

        # constructors + parent
        mol = Molecule(sys)
        prot = Protein(sys)

        # molecule_by_idx
        @test_throws KeyError molecule_by_idx(sys, -1)
        @test molecule_by_idx(sys, mol.idx) isa Molecule{T}
        @test molecule_by_idx(sys, mol.idx) == mol
        @test molecule_by_idx(sys, prot.idx) isa Molecule{T}
        @test molecule_by_idx(sys, prot.idx) == prot

        # molecules
        mt = molecules(sys)
        @test mt isa MoleculeTable{T}
        @test length(mt) == 2

        # nmolecules
        @test nmolecules(sys) isa Int
        @test nmolecules(sys) == 2

        # delete!
        frag = Fragment(Chain(mol), 1)
        Bond(Atom(frag, 1, Elements.H), Atom(frag, 2, Elements.C), BondOrder.Single)
        @test natoms(sys) == 2
        @test nbonds(sys) == 1
        @test nmolecules(sys) == 2
        @test nchains(sys) == 1
        @test nfragments(sys) == 1

        @test delete!(mol; keep_atoms = true) === nothing
        @test natoms(sys) == 2
        @test nbonds(sys) == 1
        @test nmolecules(sys) == 1
        @test nchains(sys) == 0
        @test nfragments(sys) == 0
        @test_throws KeyError mol.idx
        @test all(isnothing, atoms(sys).molecule_idx)

        frag = Fragment(Chain(prot), 1)
        Bond(Atom(frag, 1, Elements.H), Atom(frag, 2, Elements.C), BondOrder.Single)
        @test natoms(sys) == 4
        @test nbonds(sys) == 2
        @test nmolecules(sys) == 1
        @test nchains(sys) == 1
        @test nfragments(sys) == 1

        @test delete!(prot; keep_atoms = false) === nothing
        @test natoms(sys) == 2
        @test nbonds(sys) == 1
        @test nmolecules(sys) == 0
        @test nchains(sys) == 0
        @test nfragments(sys) == 0
        @test_throws KeyError prot.idx
    end
end

@testitem "Molecule/None" begin
    for T in [Float32, Float64]
        sys = System{T}()

        # constructors + parent
        mol = Molecule(sys)
        @test mol isa Molecule{T}
        @test parent(mol) === sys
        @test parent_system(mol) === sys

        if T == Float32
            mol_ds = Molecule()
            @test parent(mol_ds) === default_system()
            @test parent_system(mol_ds) === default_system()

            Molecule(; name = "something", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        mol2 = Molecule(sys; name = "something", properties = Properties(:a => 1), flags = Flags([:A, :B]))

        # getproperty
        @test mol.idx isa Int
        @test mol.name isa String
        @test mol.name == ""

        @test mol.properties isa Properties
        @test mol.properties == Properties()
        @test mol.flags isa Flags
        @test mol.flags == Flags()
        @test mol.variant isa MoleculeVariantType
        @test mol.variant == MoleculeVariant.None

        @test mol._sys isa System{T}
        @test mol._idx isa Int

        @test mol2.name == "something"
        @test mol2.properties == Properties(:a => 1)
        @test mol2.flags == Flags([:A, :B])
        @test mol2.variant == MoleculeVariant.None

        # setproperty!
        mol.name = "something else"
        @test mol.name == "something else"

        mol.properties = Properties(:first => "v1", :second => 99)
        @test length(mol.properties) == 2
        @test mol.properties[:first] == "v1"
        @test mol.properties[:second] == 99
        mol.flags = Flags([:C])
        @test length(mol.flags) == 1
        @test :C in mol.flags

        mol3 = Molecule(sys)
        mol3.variant = MoleculeVariant.Protein
        @test mol3.variant === MoleculeVariant.Protein

        # molecule_by_idx
        @test_throws KeyError molecule_by_idx(sys, -1)
        @test molecule_by_idx(sys, mol.idx) isa Molecule{T}
        @test molecule_by_idx(sys, mol.idx) == mol

        # molecules
        mt = molecules(sys; variant = MoleculeVariant.None)
        @test mt isa MoleculeTable{T}
        @test length(mt) == 2

        # nmolecules + push!
        @test nmolecules(sys; variant = MoleculeVariant.None) isa Int
        @test nmolecules(sys; variant = MoleculeVariant.None) == 2

        @test push!(sys, mol) === sys
        @test nmolecules(sys; variant = MoleculeVariant.None) == 3

        lmol = last(molecules(sys; variant = MoleculeVariant.None))
        @test lmol.idx != mol.idx
        @test lmol.name == mol.name
        @test lmol.properties == mol.properties
        @test lmol.flags == mol.flags

        # molecule atoms
        @test length(atoms(mol)) == 0
        @test atoms(mol) == atoms(sys, molecule_idx = mol.idx)
        @test natoms(mol) == 0
        @test natoms(mol) == natoms(sys, molecule_idx = mol.idx)

        @test Atom(mol, 1, Elements.H).molecule_idx == mol.idx
        @test length(atoms(mol)) == 1
        @test atoms(mol) == atoms(sys, molecule_idx = mol.idx)
        @test natoms(mol) == 1
        @test natoms(mol) == natoms(sys, molecule_idx = mol.idx)

        for atom in atoms(mol)
            @test parent_molecule(atom) === mol
            @test parent_protein(atom) === nothing
        end
        @test parent_molecule(Atom(mol, 2, Elements.C)) === mol

        # molecule bonds
        @test length(bonds(mol)) == 0
        @test bonds(mol) == bonds(sys, molecule_idx = mol.idx)
        @test nbonds(mol) == 0
        @test nbonds(mol) == nbonds(sys, molecule_idx = mol.idx)

        Bond(mol, Atom(mol, 1, Elements.H).idx, Atom(mol, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(mol)) == 1
        @test bonds(mol) == bonds(sys, molecule_idx = mol.idx)
        @test nbonds(mol) == 1
        @test nbonds(mol) == nbonds(sys, molecule_idx = mol.idx)
    end
end

@testitem "Molecule/Protein" begin
    for T in [Float32, Float64]
        sys = System{T}()

        # constructors + parent
        prot = Protein(sys)
        @test prot isa Molecule{T}
        @test parent(prot) === sys
        @test parent_system(prot) === sys

        if T == Float32
            prot_ds = Protein()
            @test parent(prot_ds) === default_system()
            @test parent_system(prot_ds) === default_system()

            Protein(; name = "something", properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        prot2 = Protein(sys; name = "something", properties = Properties(:a => 1), flags = Flags([:A, :B]))

        # getproperty
        @test prot.idx isa Int
        @test prot.name isa String
        @test prot.name == ""

        @test prot.properties isa Properties
        @test prot.properties == Properties()
        @test prot.flags isa Flags
        @test prot.flags == Flags()
        @test prot.variant isa MoleculeVariantType
        @test prot.variant === MoleculeVariant.Protein

        @test prot._sys isa System{T}
        @test prot._idx isa Int

        @test prot2.name == "something"
        @test prot2.properties == Properties(:a => 1)
        @test prot2.flags == Flags([:A, :B])
        @test prot2.variant === MoleculeVariant.Protein

        # setproperty!
        prot.name = "something else"
        @test prot.name == "something else"

        prot.properties = Properties(:first => "v1", :second => 99)
        @test length(prot.properties) == 2
        @test prot.properties[:first] == "v1"
        @test prot.properties[:second] == 99
        prot.flags = Flags([:C])
        @test length(prot.flags) == 1
        @test :C in prot.flags

        prot3 = Protein(sys)
        prot3.variant = MoleculeVariant.None
        @test prot3.variant === MoleculeVariant.None

        # protein_by_idx
        @test_throws KeyError protein_by_idx(sys, -1)
        @test protein_by_idx(sys, prot.idx) isa Molecule{T}
        @test protein_by_idx(sys, prot.idx) == prot

        # proteins
        pt = proteins(sys)
        @test pt isa MoleculeTable{T}
        @test length(pt) == 2

        # nproteins + push!
        @test nproteins(sys) isa Int
        @test nproteins(sys) == 2

        @test push!(sys, prot) === sys
        @test nproteins(sys) == 3

        lprot = last(proteins(sys))
        @test lprot.idx != prot.idx
        @test lprot.name == prot.name
        @test lprot.properties == prot.properties
        @test lprot.flags == prot.flags

        # protein atoms
        @test length(atoms(prot)) == 0
        @test atoms(prot) == atoms(sys, molecule_idx = prot.idx)
        @test natoms(prot) == 0
        @test natoms(prot) == natoms(sys, molecule_idx = prot.idx)

        @test Atom(prot, 1, Elements.H).molecule_idx == prot.idx
        @test length(atoms(prot)) == 1
        @test atoms(prot) == atoms(sys, molecule_idx = prot.idx)
        @test natoms(prot) == 1
        @test natoms(prot) == natoms(sys, molecule_idx = prot.idx)

        for atom in atoms(prot)
            @test parent_molecule(atom) === prot
            @test parent_protein(atom) === prot
        end
        @test parent_molecule(Atom(prot, 2, Elements.C)) === prot

        # molecule bonds
        @test length(bonds(prot)) == 0
        @test bonds(prot) == bonds(sys, molecule_idx = prot.idx)
        @test nbonds(prot) == 0
        @test nbonds(prot) == nbonds(sys, molecule_idx = prot.idx)

        Bond(prot, Atom(prot, 1, Elements.H).idx, Atom(prot, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(prot)) == 1
        @test bonds(prot) == bonds(sys, molecule_idx = prot.idx)
        @test nbonds(prot) == 1
        @test nbonds(prot) == nbonds(sys, molecule_idx = prot.idx)
    end
end
