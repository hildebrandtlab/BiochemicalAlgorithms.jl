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
        n1 = Nucleotide(chain, 2;
            name = "my nucleotide"
        )
        r1 = Residue(chain, 3;
            name = "my residue"
        )

        ft = fragments(sys)

        # AutoHashEquals and identity
        ft2 = fragments(sys)
        @test ft == ft2
        @test isequal(ft, ft2)
        @test hash(ft) == hash(ft2)
        @test ft !== ft2

        nt = nucleotides(sys)
        nt2 = nucleotides(sys)
        @test nt == nt2
        @test isequal(nt, nt2)
        @test hash(nt) == hash(nt2)
        @test nt !== nt2

        rt = residues(sys)
        rt2 = residues(sys)
        @test rt == rt2
        @test isequal(rt, rt2)
        @test hash(rt) == hash(rt2)
        @test rt !== rt2

        # Tables.jl interface
        @test Tables.istable(typeof(ft))
        @test Tables.columnaccess(typeof(ft))
        @test Tables.schema(ft) isa Tables.Schema
        @test !isnothing(Tables.columns(ft))
        @test !isnothing(Tables.rows(ft))

        # AbstractArray interface
        @test size(ft) == (3, 3)
        @test length(ft) == 3
        @test eltype(ft) == Fragment{T}
        @test keys(ft) == [1, 2, 3]

        # getproperty
        @test ft._sys === sys
        @test ft._idx == [f1.idx, n1.idx, r1.idx]

        @test ft.idx isa AbstractVector{Int}
        @test ft.idx == [f1.idx, n1.idx, r1.idx]
        @test ft.number isa AbstractVector{Int}
        @test ft.number == [f1.number, n1.number, r1.number]
        @test ft.name isa AbstractVector{String}
        @test ft.name == [f1.name, n1.name, r1.name]

        @test ft.properties isa AbstractVector{Properties}
        @test ft.properties == [f1.properties, n1.properties, r1.properties]
        @test ft.flags isa AbstractVector{Flags}
        @test ft.flags == [f1.flags, n1.flags, r1.flags]
        @test ft.variant isa AbstractVector{FragmentVariantType}
        @test ft.variant == [f1.variant, n1.variant, r1.variant]
        @test ft.molecule_idx isa AbstractVector{Int}
        @test ft.molecule_idx == [f1.molecule_idx, n1.molecule_idx, r1.molecule_idx]
        @test ft.chain_idx isa AbstractVector{Int}
        @test ft.chain_idx == [f1.chain_idx, n1.chain_idx, r1.chain_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(ft, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :idx) == Tables.getcolumn(ft, 1) == [f1.idx, n1.idx, r1.idx]
        @test Tables.getcolumn(ft, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :number) == Tables.getcolumn(ft, 2) == [f1.number, n1.number, r1.number]
        @test Tables.getcolumn(ft, :name) isa AbstractVector{String}
        @test Tables.getcolumn(ft, 3) isa AbstractVector{String}
        @test Tables.getcolumn(ft, :name) == Tables.getcolumn(ft, 3) == [f1.name, n1.name, r1.name]

        @test Tables.getcolumn(ft, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(ft, :properties) == [f1.properties, n1.properties, r1.properties]
        @test Tables.getcolumn(ft, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(ft, :flags) == [f1.flags, n1.flags, r1.flags]
        @test Tables.getcolumn(ft, :variant) isa AbstractVector{FragmentVariantType}
        @test Tables.getcolumn(ft, :variant) == [f1.variant, n1.variant, r1.variant]
        @test Tables.getcolumn(ft, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :molecule_idx) == [f1.molecule_idx, n1.molecule_idx, r1.molecule_idx]
        @test Tables.getcolumn(ft, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(ft, :chain_idx) == [f1.chain_idx, n1.chain_idx, r1.chain_idx]

        # setproperty!
        @test_throws ErrorException ft.idx = [999, 998, 997]
        @test_throws ErrorException ft.number = [996, 995, 994]
        @test_throws ErrorException ft.name = ["some", "other", "names"]

        @test_throws ErrorException ft.properties = [Properties(), Properties(:fourth => 995), Properties()]
        @test_throws ErrorException ft.flags = [Flags(), Flags([:fifth]), Flags()]
        @test_throws ErrorException ft.variant = [FragmentVariant.Residue, FragmentVariant.Residue, FragmentVariant.Residue]
        @test_throws ErrorException ft.molecule_idx = [993, 992, 991]
        @test_throws ErrorException ft.chain_idx = [990, 989, 988]

        # getindex
        @test ft[1] === f1
        @test ft[2] === n1
        @test ft[3] === r1
        @test_throws BoundsError ft[0]
        @test_throws BoundsError ft[4]

        # filter
        @test filter(_ -> true, ft) == ft
        @test only(filter(f -> f.idx == f1.idx, ft)) === f1

        # collect
        fv = collect(ft)
        @test fv isa Vector{Fragment{T}}
        @test length(fv) == 3
    end
end

@testitem "FragmentTable/None" begin
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

        # add decoy fragments
        Nucleotide(chain, 3)
        Residue(chain, 4)

        ft = fragments(sys; variant = FragmentVariant.None)

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
        @test all(frag -> frag.variant === FragmentVariant.None, ft)

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
        @test ft.variant isa AbstractVector{FragmentVariantType}
        @test ft.variant == [f1.variant, f2.variant]
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
        @test Tables.getcolumn(ft, :variant) isa AbstractVector{FragmentVariantType}
        @test Tables.getcolumn(ft, :variant) == [f1.variant, f2.variant]
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
        @test_throws ErrorException ft.variant = [FragmentVariant.Residue, FragmentVariant.Residue]
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

@testitem "FragmentTable/Nucleotide" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol)

        n1 = Nucleotide(chain, 1;
            name = "my nucleotide",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        n2 = Nucleotide(chain, 2)

        # add decoy fragments
        Fragment(chain, 3)
        Residue(chain, 4)

        nt = nucleotides(sys)

        # Tables.jl interface
        @test Tables.istable(typeof(nt))
        @test Tables.columnaccess(typeof(nt))
        @test Tables.schema(nt) isa Tables.Schema
        @test !isnothing(Tables.columns(nt))
        @test !isnothing(Tables.rows(nt))

        # AbstractArray interface
        @test size(nt) == (2, 3)
        @test length(nt) == 2
        @test eltype(nt) == Fragment{T}
        @test keys(nt) == [1, 2]
        @test all(isnucleotide, nt)

        # getproperty
        @test nt._sys === sys
        @test nt._idx == [n1.idx, n2.idx]

        @test nt.idx isa AbstractVector{Int}
        @test nt.idx == [n1.idx, n2.idx]
        @test nt.number isa AbstractVector{Int}
        @test nt.number == [n1.number, n2.number]
        @test nt.name isa AbstractVector{String}
        @test nt.name == [n1.name, n2.name]

        @test nt.properties isa AbstractVector{Properties}
        @test nt.properties == [n1.properties, n2.properties]
        @test nt.flags isa AbstractVector{Flags}
        @test nt.flags == [n1.flags, n2.flags]
        @test nt.variant isa AbstractVector{FragmentVariantType}
        @test nt.variant == [n1.variant, n2.variant]
        @test nt.molecule_idx isa AbstractVector{Int}
        @test nt.molecule_idx == [n1.molecule_idx, n2.molecule_idx]
        @test nt.chain_idx isa AbstractVector{Int}
        @test nt.chain_idx == [n1.chain_idx, n2.chain_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(nt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, :idx) == Tables.getcolumn(nt, 1) == [n1.idx, n2.idx]
        @test Tables.getcolumn(nt, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, :number) == Tables.getcolumn(nt, 2) == [n1.number, n2.number]
        @test Tables.getcolumn(nt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(nt, 3) isa AbstractVector{String}
        @test Tables.getcolumn(nt, :name) == Tables.getcolumn(nt, 3) == [n1.name, n2.name]

        @test Tables.getcolumn(nt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(nt, :properties) == [n1.properties, n2.properties]
        @test Tables.getcolumn(nt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(nt, :flags) == [n1.flags, n2.flags]
        @test Tables.getcolumn(nt, :variant) isa AbstractVector{FragmentVariantType}
        @test Tables.getcolumn(nt, :variant) == [n1.variant, n2.variant]
        @test Tables.getcolumn(nt, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, :molecule_idx) == [n1.molecule_idx, n2.molecule_idx]
        @test Tables.getcolumn(nt, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(nt, :chain_idx) == [n1.chain_idx, n2.chain_idx]

        # setproperty!
        @test_throws ErrorException nt.idx = [999, 998]
        @test_throws ErrorException nt.number = [997, 996]
        @test_throws ErrorException nt.name = ["some other", "names"]

        @test_throws ErrorException nt.properties = [Properties(), Properties(:fourth => 995)]
        @test_throws ErrorException nt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException nt.variant = [FragmentVariant.Residue, FragmentVariant.Residue]
        @test_throws ErrorException nt.molecule_idx = [994, 993]
        @test_throws ErrorException nt.chain_idx = [992, 991]

        # getindex
        @test nt[1] === n1
        @test nt[2] === n2
        @test_throws BoundsError nt[0]
        @test_throws BoundsError nt[3]

        # filter
        @test filter(_ -> true, nt) == nt
        @test only(filter(f -> f.idx == n1.idx, nt)) === n1

        # collect
        fv = collect(nt)
        @test fv isa Vector{Fragment{T}}
        @test length(fv) == 2
    end
end

@testitem "FragmentTable/Residue" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol)

        r1 = Residue(chain, 1;
            name = "my residue",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        r2 = Residue(chain, 2)

        # add decoy fragments
        Fragment(chain, 3)
        Nucleotide(chain, 4)

        rt = residues(sys)

        # Tables.jl interface
        @test Tables.istable(typeof(rt))
        @test Tables.columnaccess(typeof(rt))
        @test Tables.schema(rt) isa Tables.Schema
        @test !isnothing(Tables.columns(rt))
        @test !isnothing(Tables.rows(rt))

        # AbstractArray interface
        @test size(rt) == (2, 3)
        @test length(rt) == 2
        @test eltype(rt) == Fragment{T}
        @test keys(rt) == [1, 2]
        @test all(isresidue, rt)

        # getproperty
        @test rt._sys === sys
        @test rt._idx == [r1.idx, r2.idx]

        @test rt.idx isa AbstractVector{Int}
        @test rt.idx == [r1.idx, r2.idx]
        @test rt.number isa AbstractVector{Int}
        @test rt.number == [r1.number, r2.number]
        @test rt.name isa AbstractVector{String}
        @test rt.name == [r1.name, r2.name]

        @test rt.properties isa AbstractVector{Properties}
        @test rt.properties == [r1.properties, r2.properties]
        @test rt.flags isa AbstractVector{Flags}
        @test rt.flags == [r1.flags, r2.flags]
        @test rt.variant isa AbstractVector{FragmentVariantType}
        @test rt.variant == [r1.variant, r2.variant]
        @test rt.molecule_idx isa AbstractVector{Int}
        @test rt.molecule_idx == [r1.molecule_idx, r2.molecule_idx]
        @test rt.chain_idx isa AbstractVector{Int}
        @test rt.chain_idx == [r1.chain_idx, r2.chain_idx]

        # Tables.getcolumn
        @test Tables.getcolumn(rt, :idx) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, 1) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :idx) == Tables.getcolumn(rt, 1) == [r1.idx, r2.idx]
        @test Tables.getcolumn(rt, :number) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, 2) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :number) == Tables.getcolumn(rt, 2) == [r1.number, r2.number]
        @test Tables.getcolumn(rt, :name) isa AbstractVector{String}
        @test Tables.getcolumn(rt, 3) isa AbstractVector{String}
        @test Tables.getcolumn(rt, :name) == Tables.getcolumn(rt, 3) == [r1.name, r2.name]

        @test Tables.getcolumn(rt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(rt, :properties) == [r1.properties, r2.properties]
        @test Tables.getcolumn(rt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(rt, :flags) == [r1.flags, r2.flags]
        @test Tables.getcolumn(rt, :variant) isa AbstractVector{FragmentVariantType}
        @test Tables.getcolumn(rt, :variant) == [r1.variant, r2.variant]
        @test Tables.getcolumn(rt, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :molecule_idx) == [r1.molecule_idx, r2.molecule_idx]
        @test Tables.getcolumn(rt, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :chain_idx) == [r1.chain_idx, r2.chain_idx]

        # setproperty!
        @test_throws ErrorException rt.idx = [999, 998]
        @test_throws ErrorException rt.number = [997, 996]
        @test_throws ErrorException rt.name = ["some other", "names"]

        @test_throws ErrorException rt.properties = [Properties(), Properties(:fourth => 995)]
        @test_throws ErrorException rt.flags = [Flags(), Flags([:fifth])]
        @test_throws ErrorException rt.variant = [FragmentVariant.None, FragmentVariant.None]
        @test_throws ErrorException rt.molecule_idx = [994, 993]
        @test_throws ErrorException rt.chain_idx = [992, 991]

        # getindex
        @test rt[1] === r1
        @test rt[2] === r2
        @test_throws BoundsError rt[0]
        @test_throws BoundsError rt[3]

        # filter
        @test filter(_ -> true, rt) == rt
        @test only(filter(f -> f.idx == r1.idx, rt)) === r1

        # collect
        fv = collect(rt)
        @test fv isa Vector{Fragment{T}}
        @test length(fv) == 2
    end
end

@testitem "Fragment" begin
    for T in [Float32, Float64]
        sys = System{T}()
        chain = Chain(Molecule(sys))

        # constructors + parent
        frag = Fragment(chain, 1)
        nuc = Nucleotide(chain, 2)
        res = Residue(chain, 3)

        # fragment_by_idx
        @test_throws KeyError fragment_by_idx(sys, -1)
        @test fragment_by_idx(sys, frag.idx) isa Fragment{T}
        @test fragment_by_idx(sys, frag.idx) == frag
        @test fragment_by_idx(sys, nuc.idx) isa Fragment{T}
        @test fragment_by_idx(sys, nuc.idx) == nuc
        @test fragment_by_idx(sys, res.idx) isa Fragment{T}
        @test fragment_by_idx(sys, res.idx) == res

        # fragments
        ft = fragments(sys)
        @test ft isa FragmentTable{T}
        @test length(ft) == 3
        @test length(fragments(sys)) == 3

        # nfragments
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 3
    end
end

@testitem "Fragment/None" begin
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
        @test length(frag._row) == 3

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
        @test frag.variant isa FragmentVariantType
        @test frag.variant === FragmentVariant.None
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
        @test frag2.variant == FragmentVariant.None
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

        frag3 = Fragment(Chain(Molecule(sys)), 1)
        frag3.variant = FragmentVariant.Residue
        @test frag3.variant === FragmentVariant.Residue
        frag3.molecule_idx = 999
        @test frag3.molecule_idx == 999
        frag3.chain_idx = 998
        @test frag3.chain_idx == 998

        # fragment_by_idx
        @test_throws KeyError fragment_by_idx(sys, -1)
        @test fragment_by_idx(sys, frag.idx) isa Fragment{T}
        @test fragment_by_idx(sys, frag.idx) == frag

        # fragments
        ft = fragments(sys; variant = FragmentVariant.None)
        @test ft isa FragmentTable{T}
        @test length(ft) == 2
        @test length(fragments(sys; variant = FragmentVariant.None)) == 2
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = -1)) == 0
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = mol2.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = nothing)) == 2
        @test length(fragments(sys; variant = FragmentVariant.None, chain_idx = -1)) == 0
        @test length(fragments(sys; variant = FragmentVariant.None, chain_idx = chain.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, chain_idx = chain2.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, chain_idx = nothing)) == 2
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = -1, chain_idx = chain.idx)) == 0
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = -1)) == 0
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = chain.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = nothing)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = chain.idx)) == 1
        @test length(fragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = nothing)) == 2

        # nfragments + push!
        @test nfragments(sys; variant = FragmentVariant.None) isa Int
        @test nfragments(sys; variant = FragmentVariant.None) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol2.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing) == 2
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = chain.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = chain2.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = nothing) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = chain.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = nothing) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = chain.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = nothing) == 2

        @test push!(chain, frag) === chain
        @test nfragments(sys; variant = FragmentVariant.None) isa Int
        @test nfragments(sys; variant = FragmentVariant.None) == 3
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol2.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing) == 3
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = chain.idx) == 2
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = chain2.idx) == 1
        @test nfragments(sys; variant = FragmentVariant.None, chain_idx = nothing) == 3
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = chain.idx) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol.idx, chain_idx = nothing) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = chain.idx) == 2
        @test nfragments(sys; variant = FragmentVariant.None, molecule_idx = nothing, chain_idx = nothing) == 3

        lfrag = last(fragments(chain; variant = FragmentVariant.None))
        @test lfrag.idx != frag.idx
        @test lfrag.number == frag.number
        @test lfrag.name == frag.name
        @test lfrag.properties == frag.properties
        @test lfrag.flags == frag.flags
        @test lfrag.variant == frag.variant
        @test lfrag.molecule_idx == parent_molecule(chain).idx
        @test lfrag.chain_idx == chain.idx

        # chain/molecule fragments
        mol3 = Molecule(sys)
        @test size(fragments(mol3; variant = FragmentVariant.None), 1) == 0
        @test fragments(mol3) == fragments(sys; variant = FragmentVariant.None, molecule_idx = mol3.idx)
        @test nfragments(mol3) == 0
        @test nfragments(mol3) == nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol3.idx)

        chain3 = Chain(mol3)
        @test size(fragments(chain3; variant = FragmentVariant.None), 1) == 0
        @test fragments(chain3; variant = FragmentVariant.None) == fragments(sys; variant = FragmentVariant.None, chain_idx = chain3.idx)
        @test nfragments(chain3; variant = FragmentVariant.None) == 0
        @test nfragments(chain3; variant = FragmentVariant.None) == nfragments(sys; variant = FragmentVariant.None, chain_idx = chain3.idx)

        Fragment(chain3, 1)
        @test size(fragments(mol3; variant = FragmentVariant.None), 1) == 1
        @test fragments(mol3; variant = FragmentVariant.None) == fragments(sys; variant = FragmentVariant.None, molecule_idx = mol3.idx)
        @test nfragments(mol3; variant = FragmentVariant.None) == 1
        @test nfragments(mol3; variant = FragmentVariant.None) == nfragments(sys; variant = FragmentVariant.None, molecule_idx = mol3.idx)

        @test size(fragments(chain3; variant = FragmentVariant.None), 1) == 1
        @test fragments(chain3; variant = FragmentVariant.None) == fragments(sys; variant = FragmentVariant.None, chain_idx = chain3.idx)
        @test nfragments(chain3; variant = FragmentVariant.None) == 1
        @test nfragments(chain3; variant = FragmentVariant.None) == nfragments(sys; variant = FragmentVariant.None, chain_idx = chain3.idx)

        # fragment atoms
        @test length(atoms(frag)) == 0
        @test atoms(frag) == atoms(sys; fragment_idx = frag.idx)
        @test natoms(frag) == 0
        @test natoms(frag) == natoms(sys; fragment_idx = frag.idx)

        @test Atom(frag, 1, Elements.H).fragment_idx == frag.idx
        @test length(atoms(frag)) == 1
        @test atoms(frag) == atoms(sys; fragment_idx = frag.idx)
        @test natoms(frag) == 1
        @test natoms(frag) == natoms(sys; fragment_idx = frag.idx)

        for atom in atoms(frag)
            @test parent_fragment(atom) === frag
            @test parent_nucleotide(atom) === nothing
            @test parent_residue(atom) === nothing
        end

        # fragment bonds
        @test length(bonds(frag)) == 0
        @test bonds(frag) == bonds(sys; fragment_idx = frag.idx)
        @test nbonds(frag) == 0
        @test nbonds(frag) == nbonds(sys; fragment_idx = frag.idx)

        Bond(frag, Atom(frag, 1, Elements.H).idx, Atom(frag, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(frag)) == 1
        @test bonds(frag) == bonds(sys; fragment_idx = frag.idx)
        @test nbonds(frag) == 1
        @test nbonds(frag) == nbonds(sys; fragment_idx = frag.idx)
    end
end

@testitem "Fragment/Nucleotide" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        nuc = Nucleotide(chain, 1)
        @test nuc isa Fragment{T}
        @test parent(nuc) === sys
        @test parent_system(nuc) === sys
        @test parent_molecule(nuc) === mol
        @test parent_chain(nuc) === chain

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            nuc_ds = Nucleotide(chain_ds, 1)
            parent(nuc_ds) === default_system()
            parent_system(nuc_ds) === default_system()

            Nucleotide(chain_ds, 1;
                name = "something",
                properties = Properties(:a => "b"),
                flags = Flags([:A])
            )
        end

        nuc2 = Nucleotide(chain2, 1;
            name = "something",
            properties = Properties(:a => 1),
            flags = Flags([:A, :B])
        )

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(nuc._row) == 3

        # getproperty
        @test nuc.idx isa Int
        @test nuc.number isa Int
        @test nuc.number == 1
        @test nuc.name isa String
        @test nuc.name == ""

        @test nuc.properties isa Properties
        @test nuc.properties == Properties()
        @test nuc.flags isa Flags
        @test nuc.flags == Flags()
        @test nuc.variant isa FragmentVariantType
        @test nuc.variant === FragmentVariant.Nucleotide
        @test nuc.molecule_idx isa Int
        @test nuc.molecule_idx == mol.idx
        @test nuc.chain_idx isa Int
        @test nuc.chain_idx == chain.idx

        @test nuc._sys isa System{T}
        @test nuc._row isa BiochemicalAlgorithms._FragmentTableRow
        
        @test nuc2.number == 1
        @test nuc2.name == "something"
        @test nuc2.properties == Properties(:a => 1)
        @test nuc2.flags == Flags([:A, :B])
        @test nuc2.variant === FragmentVariant.Nucleotide
        @test nuc2.molecule_idx == mol2.idx
        @test nuc2.chain_idx == chain2.idx

        # setproperty!
        nuc.number = 0
        @test nuc.number == 0
        nuc.name = "something else"
        @test nuc.name == "something else"

        nuc.properties = Properties(:first => "v1", :second => 99)
        @test length(nuc.properties) == 2
        @test nuc.properties[:first] == "v1"
        @test nuc.properties[:second] == 99
        nuc.flags = Flags([:C])
        @test length(nuc.flags) == 1
        @test :C in nuc.flags

        nuc3 = Nucleotide(Chain(Molecule()), 1)
        nuc3.variant = FragmentVariant.None
        @test nuc3.variant === FragmentVariant.None
        nuc3.molecule_idx = 999
        @test nuc3.molecule_idx == 999
        nuc3.chain_idx = 998
        @test nuc3.chain_idx == 998

        # nucleotide_by_idx
        @test_throws KeyError nucleotide_by_idx(sys, -1)
        @test nucleotide_by_idx(sys, nuc.idx) isa Fragment{T}
        @test nucleotide_by_idx(sys, nuc.idx) == nuc

        # nucleotides
        nt = nucleotides(sys)
        @test nt isa FragmentTable{T}
        @test length(nt) == 2
        @test length(nucleotides(sys)) == 2
        @test length(nucleotides(sys, molecule_idx = -1)) == 0
        @test length(nucleotides(sys, molecule_idx = mol.idx)) == 1
        @test length(nucleotides(sys, molecule_idx = mol2.idx)) == 1
        @test length(nucleotides(sys, molecule_idx = nothing)) == 2
        @test length(nucleotides(sys, chain_idx = -1)) == 0
        @test length(nucleotides(sys, chain_idx = chain.idx)) == 1
        @test length(nucleotides(sys, chain_idx = chain2.idx)) == 1
        @test length(nucleotides(sys, chain_idx = nothing)) == 2
        @test length(nucleotides(sys, molecule_idx = -1, chain_idx = chain.idx)) == 0
        @test length(nucleotides(sys, molecule_idx = mol.idx, chain_idx = -1)) == 0
        @test length(nucleotides(sys, molecule_idx = mol.idx, chain_idx = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_idx = mol.idx, chain_idx = nothing)) == 1
        @test length(nucleotides(sys, molecule_idx = nothing, chain_idx = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_idx = nothing, chain_idx = nothing)) == 2

        # nnucleotides + push!
        @test nnucleotides(sys) isa Int
        @test nnucleotides(sys) == 2
        @test nnucleotides(sys, molecule_idx = -1) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx) == 1
        @test nnucleotides(sys, molecule_idx = mol2.idx) == 1
        @test nnucleotides(sys, molecule_idx = nothing) == 2
        @test nnucleotides(sys, chain_idx = -1) == 0
        @test nnucleotides(sys, chain_idx = chain.idx) == 1
        @test nnucleotides(sys, chain_idx = chain2.idx) == 1
        @test nnucleotides(sys, chain_idx = nothing) == 2
        @test nnucleotides(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 1
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = nothing) == 1
        @test nnucleotides(sys, molecule_idx = nothing, chain_idx = chain.idx) == 1
        @test nnucleotides(sys, molecule_idx = nothing, chain_idx = nothing) == 2

        @test push!(chain, nuc) === chain
        @test nnucleotides(sys) isa Int
        @test nnucleotides(sys) == 3
        @test nnucleotides(sys, molecule_idx = -1) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx) == 2
        @test nnucleotides(sys, molecule_idx = mol2.idx) == 1
        @test nnucleotides(sys, molecule_idx = nothing) == 3
        @test nnucleotides(sys, chain_idx = -1) == 0
        @test nnucleotides(sys, chain_idx = chain.idx) == 2
        @test nnucleotides(sys, chain_idx = chain2.idx) == 1
        @test nnucleotides(sys, chain_idx = nothing) == 3
        @test nnucleotides(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 2
        @test nnucleotides(sys, molecule_idx = mol.idx, chain_idx = nothing) == 2
        @test nnucleotides(sys, molecule_idx = nothing, chain_idx = chain.idx) == 2
        @test nnucleotides(sys, molecule_idx = nothing, chain_idx = nothing) == 3

        lnuc = last(nucleotides(chain))
        @test lnuc.idx != nuc.idx
        @test lnuc.number == nuc.number
        @test lnuc.name == nuc.name
        @test lnuc.properties == nuc.properties
        @test lnuc.flags == nuc.flags
        @test lnuc.variant == nuc.variant
        @test lnuc.molecule_idx == parent_molecule(chain).idx
        @test lnuc.chain_idx == chain.idx

        # chain/molecule nucleotides
        mol3 = Molecule(sys)
        @test size(nucleotides(mol3), 1) == 0
        @test nucleotides(mol3) == nucleotides(sys, molecule_idx = mol3.idx)
        @test nnucleotides(mol3) == 0
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_idx = mol3.idx)

        chain3 = Chain(mol3)
        @test size(nucleotides(chain3), 1) == 0
        @test nucleotides(chain3) == nucleotides(sys, chain_idx = chain3.idx)
        @test nnucleotides(chain3) == 0
        @test nnucleotides(chain3) == nnucleotides(sys, chain_idx = chain3.idx)

        Nucleotide(chain3, 1)
        @test size(nucleotides(mol3), 1) == 1
        @test nucleotides(mol3) == nucleotides(sys, molecule_idx = mol3.idx)
        @test nnucleotides(mol3) == 1
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_idx = mol3.idx)

        @test size(nucleotides(chain3), 1) == 1
        @test nucleotides(chain3) == nucleotides(sys, chain_idx = chain3.idx)
        @test nnucleotides(chain3) == 1
        @test nnucleotides(chain3) == nnucleotides(sys, chain_idx = chain3.idx)

        # nucleotide atoms
        @test length(atoms(nuc)) == 0
        @test atoms(nuc) == atoms(sys, fragment_idx = nuc.idx)
        @test natoms(nuc) == 0
        @test natoms(nuc) == natoms(sys, fragment_idx = nuc.idx)

        @test Atom(nuc, 1, Elements.H).fragment_idx == nuc.idx
        @test length(atoms(nuc)) == 1
        @test atoms(nuc) == atoms(sys, fragment_idx = nuc.idx)
        @test natoms(nuc) == 1
        @test natoms(nuc) == natoms(sys, fragment_idx = nuc.idx)

        for atom in atoms(nuc)
            @test parent_fragment(atom) === nuc
            @test parent_nucleotide(atom) === nuc
            @test parent_residue(atom) === nothing
        end

        # nucleotide bonds
        @test length(bonds(nuc)) == 0
        @test bonds(nuc) == bonds(sys, fragment_idx = nuc.idx)
        @test nbonds(nuc) == 0
        @test nbonds(nuc) == nbonds(sys, fragment_idx = nuc.idx)

        Bond(nuc, Atom(nuc, 1, Elements.H).idx, Atom(nuc, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(nuc)) == 1
        @test bonds(nuc) == bonds(sys, fragment_idx = nuc.idx)
        @test nbonds(nuc) == 1
        @test nbonds(nuc) == nbonds(sys, fragment_idx = nuc.idx)
    end
end

@testitem "Fragment/Residue" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        res = Residue(chain, 1)
        @test res isa Fragment{T}
        @test parent(res) === sys
        @test parent_system(res) === sys
        @test parent_molecule(res) === mol
        @test parent_chain(res) === chain

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            res_ds = Residue(chain_ds, 1)
            parent(res_ds) === default_system()
            parent_system(res_ds) === default_system()

            Residue(chain_ds, 1;
                name = "something",
                properties = Properties(:a => "b"),
                flags = Flags([:A])
            )
        end

        res2 = Residue(chain2, 1;
            name = "something",
            properties = Properties(:a => 1),
            flags = Flags([:A, :B])
        )

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(res._row) == 3

        # getproperty
        @test res.idx isa Int
        @test res.number isa Int
        @test res.number == 1
        @test res.name isa String
        @test res.name == ""

        @test res.properties isa Properties
        @test res.properties == Properties()
        @test res.flags isa Flags
        @test res.flags == Flags()
        @test res.variant isa FragmentVariantType
        @test res.variant === FragmentVariant.Residue
        @test res.molecule_idx isa Int
        @test res.molecule_idx == mol.idx
        @test res.chain_idx isa Int
        @test res.chain_idx == chain.idx

        @test res._sys isa System{T}
        @test res._row isa BiochemicalAlgorithms._FragmentTableRow
        
        @test res2.number == 1
        @test res2.name == "something"
        @test res2.properties == Properties(:a => 1)
        @test res2.flags == Flags([:A, :B])
        @test res2.variant === FragmentVariant.Residue
        @test res2.molecule_idx == mol2.idx
        @test res2.chain_idx == chain2.idx

        # setproperty!
        res.number = 0
        @test res.number == 0
        res.name = "something else"
        @test res.name == "something else"

        res.properties = Properties(:first => "v1", :second => 99)
        @test length(res.properties) == 2
        @test res.properties[:first] == "v1"
        @test res.properties[:second] == 99
        res.flags = Flags([:C])
        @test length(res.flags) == 1
        @test :C in res.flags

        res3 = Residue(Chain(Molecule()), 1)
        res3.variant = FragmentVariant.None
        @test res3.variant === FragmentVariant.None
        res3.molecule_idx = 999
        @test res3.molecule_idx == 999
        res3.chain_idx = 998
        @test res3.chain_idx == 998

        # residue_by_idx
        @test_throws KeyError residue_by_idx(sys, -1)
        @test residue_by_idx(sys, res.idx) isa Fragment{T}
        @test residue_by_idx(sys, res.idx) == res

        # residues
        rt = residues(sys)
        @test rt isa FragmentTable{T}
        @test length(rt) == 2
        @test length(residues(sys)) == 2
        @test length(residues(sys, molecule_idx = -1)) == 0
        @test length(residues(sys, molecule_idx = mol.idx)) == 1
        @test length(residues(sys, molecule_idx = mol2.idx)) == 1
        @test length(residues(sys, molecule_idx = nothing)) == 2
        @test length(residues(sys, chain_idx = -1)) == 0
        @test length(residues(sys, chain_idx = chain.idx)) == 1
        @test length(residues(sys, chain_idx = chain2.idx)) == 1
        @test length(residues(sys, chain_idx = nothing)) == 2
        @test length(residues(sys, molecule_idx = -1, chain_idx = chain.idx)) == 0
        @test length(residues(sys, molecule_idx = mol.idx, chain_idx = -1)) == 0
        @test length(residues(sys, molecule_idx = mol.idx, chain_idx = chain.idx)) == 1
        @test length(residues(sys, molecule_idx = mol.idx, chain_idx = nothing)) == 1
        @test length(residues(sys, molecule_idx = nothing, chain_idx = chain.idx)) == 1
        @test length(residues(sys, molecule_idx = nothing, chain_idx = nothing)) == 2

        # nresidues + push!
        @test nresidues(sys) isa Int
        @test nresidues(sys) == 2
        @test nresidues(sys, molecule_idx = -1) == 0
        @test nresidues(sys, molecule_idx = mol.idx) == 1
        @test nresidues(sys, molecule_idx = mol2.idx) == 1
        @test nresidues(sys, molecule_idx = nothing) == 2
        @test nresidues(sys, chain_idx = -1) == 0
        @test nresidues(sys, chain_idx = chain.idx) == 1
        @test nresidues(sys, chain_idx = chain2.idx) == 1
        @test nresidues(sys, chain_idx = nothing) == 2
        @test nresidues(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 1
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = nothing) == 1
        @test nresidues(sys, molecule_idx = nothing, chain_idx = chain.idx) == 1
        @test nresidues(sys, molecule_idx = nothing, chain_idx = nothing) == 2

        @test push!(chain, res) === chain
        @test nresidues(sys) isa Int
        @test nresidues(sys) == 3
        @test nresidues(sys, molecule_idx = -1) == 0
        @test nresidues(sys, molecule_idx = mol.idx) == 2
        @test nresidues(sys, molecule_idx = mol2.idx) == 1
        @test nresidues(sys, molecule_idx = nothing) == 3
        @test nresidues(sys, chain_idx = -1) == 0
        @test nresidues(sys, chain_idx = chain.idx) == 2
        @test nresidues(sys, chain_idx = chain2.idx) == 1
        @test nresidues(sys, chain_idx = nothing) == 3
        @test nresidues(sys, molecule_idx = -1, chain_idx = chain.idx) == 0
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = -1) == 0
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = chain.idx) == 2
        @test nresidues(sys, molecule_idx = mol.idx, chain_idx = nothing) == 2
        @test nresidues(sys, molecule_idx = nothing, chain_idx = chain.idx) == 2
        @test nresidues(sys, molecule_idx = nothing, chain_idx = nothing) == 3

        lres = last(residues(chain))
        @test lres.idx != res.idx
        @test lres.number == res.number
        @test lres.name == res.name
        @test lres.properties == res.properties
        @test lres.flags == res.flags
        @test lres.variant == res.variant
        @test lres.molecule_idx == parent_molecule(chain).idx
        @test lres.chain_idx == chain.idx

        # chain/molecule residues
        mol3 = Molecule(sys)
        @test size(residues(mol3), 1) == 0
        @test residues(mol3) == residues(sys, molecule_idx = mol3.idx)
        @test nresidues(mol3) == 0
        @test nresidues(mol3) == nresidues(sys, molecule_idx = mol3.idx)

        chain3 = Chain(mol3)
        @test size(residues(chain3), 1) == 0
        @test residues(chain3) == residues(sys, chain_idx = chain3.idx)
        @test nresidues(chain3) == 0
        @test nresidues(chain3) == nresidues(sys, chain_idx = chain3.idx)

        Residue(chain3, 1)
        @test size(residues(mol3), 1) == 1
        @test residues(mol3) == residues(sys, molecule_idx = mol3.idx)
        @test nresidues(mol3) == 1
        @test nresidues(mol3) == nresidues(sys, molecule_idx = mol3.idx)

        @test size(residues(chain3), 1) == 1
        @test residues(chain3) == residues(sys, chain_idx = chain3.idx)
        @test nresidues(chain3) == 1
        @test nresidues(chain3) == nresidues(sys, chain_idx = chain3.idx)

        # nucleotide atoms
        @test length(atoms(res)) == 0
        @test atoms(res) == atoms(sys, fragment_idx = res.idx)
        @test natoms(res) == 0
        @test natoms(res) == natoms(sys, fragment_idx = res.idx)

        @test Atom(res, 1, Elements.H).fragment_idx == res.idx
        @test length(atoms(res)) == 1
        @test atoms(res) == atoms(sys, fragment_idx = res.idx)
        @test natoms(res) == 1
        @test natoms(res) == natoms(sys, fragment_idx = res.idx)

        for atom in atoms(res)
            @test parent_fragment(atom) === res
            @test parent_nucleotide(atom) === nothing
            @test parent_residue(atom) === res
        end

        # nucleotide bonds
        @test length(bonds(res)) == 0
        @test bonds(res) == bonds(sys, fragment_idx = res.idx)
        @test nbonds(res) == 0
        @test nbonds(res) == nbonds(sys, fragment_idx = res.idx)

        Bond(res, Atom(res, 1, Elements.H).idx, Atom(res, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(res)) == 1
        @test bonds(res) == bonds(sys, fragment_idx = res.idx)
        @test nbonds(res) == 1
        @test nbonds(res) == nbonds(sys, fragment_idx = res.idx)
    end
end
