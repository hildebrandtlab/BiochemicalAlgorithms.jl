@testitem "NucleotideTable" begin
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
        @test eltype(nt) == Nucleotide{T}
        @test keys(nt) == [1, 2]

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
        @test_throws ErrorException nt.molecule_idx = [994, 993]
        @test_throws ErrorException nt.chain_idx = [992, 991]

        # getindex
        @test nt[1] === n1
        @test nt[2] === n2
        @test_throws BoundsError nt[0]
        @test_throws BoundsError nt[3]

        # filter
        @test filter(_ -> true, nt) == nt
        @test only(filter(n -> n.idx == n1.idx, nt)) === n1

        # collect
        nv = collect(nt)
        @test nv isa Vector{Nucleotide{T}}
        @test length(nv) == 2
    end
end

@testitem "Nucleotide" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        nuc = Nucleotide(chain, 1)
        @test nuc isa Nucleotide{T}
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
        @test length(getfield(nuc, :_row)) == 3

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
        @test nuc.molecule_idx isa Int
        @test nuc.molecule_idx == mol.idx
        @test nuc.chain_idx isa Int
        @test nuc.chain_idx == chain.idx

        @test nuc._sys isa System{T}
        @test nuc._row isa BiochemicalAlgorithms._NucleotideTableRow

        @test nuc2.number == 1
        @test nuc2.name == "something"
        @test nuc2.properties == Properties(:a => 1)
        @test nuc2.flags == Flags([:A, :B])
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

        nuc3 = Nucleotide(Chain(Molecule(System{T}())), 1)
        nuc3.molecule_idx = 999
        @test nuc3.molecule_idx == 999
        nuc3.chain_idx = 998
        @test nuc3.chain_idx == 998

        # nucleotide_by_idx
        @test_throws KeyError nucleotide_by_idx(sys, -1)
        @test nucleotide_by_idx(sys, nuc.idx) isa Nucleotide{T}
        @test nucleotide_by_idx(sys, nuc.idx) == nuc

        # nucleotides
        fv = nucleotides(sys)
        @test fv isa NucleotideTable{T}
        @test length(fv) == 2
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
        @test atoms(nuc) == atoms(sys, nucleotide_idx = nuc.idx)
        @test natoms(nuc) == 0
        @test natoms(nuc) == natoms(sys, nucleotide_idx = nuc.idx)

        @test Atom(nuc, 1, Elements.H).nucleotide_idx == nuc.idx
        @test length(atoms(nuc)) == 1
        @test atoms(nuc) == atoms(sys, nucleotide_idx = nuc.idx)
        @test natoms(nuc) == 1
        @test natoms(nuc) == natoms(sys, nucleotide_idx = nuc.idx)

        for atom in atoms(nuc)
            @test parent_nucleotide(atom) === nuc
        end

        # nucleotide bonds
        @test length(bonds(nuc)) == 0
        @test bonds(nuc) == bonds(sys, nucleotide_idx = nuc.idx)
        @test nbonds(nuc) == 0
        @test nbonds(nuc) == nbonds(sys, nucleotide_idx = nuc.idx)

        Bond(nuc, Atom(nuc, 1, Elements.H).idx, Atom(nuc, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(nuc)) == 1
        @test bonds(nuc) == bonds(sys, nucleotide_idx = nuc.idx)
        @test nbonds(nuc) == 1
        @test nbonds(nuc) == nbonds(sys, nucleotide_idx = nuc.idx)
    end
end
