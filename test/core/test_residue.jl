@testitem "ResidueTable" begin
    using Tables

    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        chain = Chain(mol)

        r1 = Residue(chain, 1, AminoAcid('A');
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        r2 = Residue(chain, 2, AminoAcid('D'))

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
        @test eltype(rt) == Residue{T}
        @test keys(rt) == [1, 2]

        # getproperty
        @test rt._sys === sys
        @test rt._idx == [r1.idx, r2.idx]

        @test rt.idx isa AbstractVector{Int}
        @test rt.idx == [r1.idx, r2.idx]
        @test rt.number isa AbstractVector{Int}
        @test rt.number == [r1.number, r2.number]
        @test rt.type isa AbstractVector{AminoAcid}
        @test rt.type == [r1.type, r2.type]

        @test rt.properties isa AbstractVector{Properties}
        @test rt.properties == [r1.properties, r2.properties]
        @test rt.flags isa AbstractVector{Flags}
        @test rt.flags == [r1.flags, r2.flags]
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
        @test Tables.getcolumn(rt, :type) isa AbstractVector{AminoAcid}
        @test Tables.getcolumn(rt, 3) isa AbstractVector{AminoAcid}
        @test Tables.getcolumn(rt, :type) == Tables.getcolumn(rt, 3) == [r1.type, r2.type]

        @test Tables.getcolumn(rt, :properties) isa AbstractVector{Properties}
        @test Tables.getcolumn(rt, :properties) == [r1.properties, r2.properties]
        @test Tables.getcolumn(rt, :flags) isa AbstractVector{Flags}
        @test Tables.getcolumn(rt, :flags) == [r1.flags, r2.flags]
        @test Tables.getcolumn(rt, :molecule_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :molecule_idx) == [r1.molecule_idx, r2.molecule_idx]
        @test Tables.getcolumn(rt, :chain_idx) isa AbstractVector{Int}
        @test Tables.getcolumn(rt, :chain_idx) == [r1.chain_idx, r2.chain_idx]

        # setproperty!
        @test_throws ErrorException rt.idx = [999, 998]
        @test_throws ErrorException rt.number = [997, 996]
        @test_throws ErrorException rt.type = [AminoAcid('D'), AminoAcid('A')]

        @test_throws ErrorException rt.properties = [Properties(), Properties(:fourth => 995)]
        @test_throws ErrorException rt.flags = [Flags(), Flags([:fifth])]
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
        @test fv isa Vector{Residue{T}}
        @test length(fv) == 2
    end
end

@testitem "Residue" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        res = Residue(chain, 1, AminoAcid('A'))
        @test res isa Residue{T}
        @test parent(res) === sys
        @test parent_system(res) === sys
        @test parent_molecule(res) === mol
        @test parent_chain(res) === chain

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            res_ds = Residue(chain_ds, 1, AminoAcid('A'))
            parent(res_ds) === default_system()
            parent_system(res_ds) === default_system()

            Residue(chain_ds, 1, AminoAcid('D'); properties = Properties(:a => "b"), flags = Flags([:A]))
        end

        res2 = Residue(chain2, 1, AminoAcid('D'); properties = Properties(:a => 1), flags = Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(res, :_row)) == 3

        # getproperty
        @test res.idx isa Int
        @test res.number isa Int
        @test res.number == 1
        @test res.type isa AminoAcid
        @test res.type == AminoAcid('A')

        @test res.properties isa Properties
        @test res.properties == Properties()
        @test res.flags isa Flags
        @test res.flags == Flags()
        @test res.molecule_idx isa Int
        @test res.molecule_idx == mol.idx
        @test res.chain_idx isa Int
        @test res.chain_idx == chain.idx

        @test res._sys isa System{T}
        @test res._row isa BiochemicalAlgorithms._ResidueTableRow

        @test res2.number == 1
        @test res2.type == AminoAcid('D')
        @test res2.properties == Properties(:a => 1)
        @test res2.flags == Flags([:A, :B])
        @test res2.molecule_idx == mol2.idx
        @test res2.chain_idx == chain2.idx

        # setproperty!
        res.number = 0
        @test res.number == 0
        res.type = AminoAcid('W')
        @test res.type == AminoAcid('W')

        res.properties = Properties(:first => "v1", :second => 99)
        @test length(res.properties) == 2
        @test res.properties[:first] == "v1"
        @test res.properties[:second] == 99
        res.flags = Flags([:C])
        @test length(res.flags) == 1
        @test :C in res.flags

        res3 = Nucleotide(Chain(Molecule(System{T}())), 1)
        res3.molecule_idx = 999
        @test res3.molecule_idx == 999
        res3.chain_idx = 998
        @test res3.chain_idx == 998

        # residue_by_idx
        @test_throws KeyError residue_by_idx(sys, -1)
        @test residue_by_idx(sys, res.idx) isa Residue{T}
        @test residue_by_idx(sys, res.idx) == res

        # residues
        fv = residues(sys)
        @test fv isa ResidueTable{T}
        @test length(fv) == 2
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
        @test lres.type == res.type
        @test lres.properties == res.properties
        @test lres.flags == res.flags
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

        Residue(chain3, 1, AminoAcid('A'))
        @test size(residues(mol3), 1) == 1
        @test residues(mol3) == residues(sys, molecule_idx = mol3.idx)
        @test nresidues(mol3) == 1
        @test nresidues(mol3) == nresidues(sys, molecule_idx = mol3.idx)

        @test size(residues(chain3), 1) == 1
        @test residues(chain3) == residues(sys, chain_idx = chain3.idx)
        @test nresidues(chain3) == 1
        @test nresidues(chain3) == nresidues(sys, chain_idx = chain3.idx)

        # residue atoms
        @test length(atoms(res)) == 0
        @test atoms(res) == atoms(sys, residue_idx = res.idx)
        @test natoms(res) == 0
        @test natoms(res) == natoms(sys, residue_idx = res.idx)

        @test Atom(res, 1, Elements.H).residue_idx == res.idx
        @test length(atoms(res)) == 1
        @test atoms(res) == atoms(sys, residue_idx = res.idx)
        @test natoms(res) == 1
        @test natoms(res) == natoms(sys, residue_idx = res.idx)

        for atom in atoms(res)
            @test parent_residue(atom) === res
        end

        # residue bonds
        @test length(bonds(res)) == 0
        @test bonds(res) == bonds(sys, residue_idx = res.idx)
        @test nbonds(res) == 0
        @test nbonds(res) == nbonds(sys, residue_idx = res.idx)

        Bond(res, Atom(res, 1, Elements.H).idx, Atom(res, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(res)) == 1
        @test bonds(res) == bonds(sys, residue_idx = res.idx)
        @test nbonds(res) == 1
        @test nbonds(res) == nbonds(sys, residue_idx = res.idx)
    end
end
