@testitem "Residue" begin
    using BiochemicalAlgorithms: _SystemResidueTuple, _residues, _atoms, _bonds
    using DataFrames

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

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            res_ds = Residue(chain_ds, 1, AminoAcid('A'))
            parent(res_ds) === default_system()
            parent_system(res_ds) === default_system()

            Residue(chain_ds, 1, AminoAcid('D'), Properties(:a => "b"), Flags([:A]))
        end

        res2 = Residue(chain2, 1, AminoAcid('D'), Properties(:a => 1), Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(res, :_row)) == 7

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

        @test res._sys isa System{T}
        @test res._row isa DataFrameRow
        
        @test_throws ErrorException res.molecule_id
        @test_throws ErrorException res.chain_id

        @test res._row.molecule_id isa Int
        @test res._row.molecule_id == mol.idx
        @test res._row.chain_id isa Int
        @test res._row.chain_id == chain.idx

        @test res2.number == 1
        @test res2.type == AminoAcid('D')
        @test res2.properties == Properties(:a => 1)
        @test res2.flags == Flags([:A, :B])
        @test res2._row.molecule_id == mol2.idx
        @test res2._row.chain_id == chain2.idx

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

        @test_throws ErrorException res.molecule_id = 0
        @test_throws ErrorException res.chain_id = 0

        # residue_by_idx
        @test_throws KeyError residue_by_idx(sys, -1)
        @test residue_by_idx(sys, res.idx) isa Residue{T}
        @test residue_by_idx(sys, res.idx) == res

        # _residues
        df = _residues(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(_SystemResidueTuple)))
        @test copy(df[1, 1:length(fieldnames(ResidueTuple))]) isa ResidueTuple
        @test size(_residues(sys), 1) == 2
        @test size(_residues(sys, molecule_id = -1), 1) == 0
        @test size(_residues(sys, molecule_id = mol.idx), 1) == 1
        @test size(_residues(sys, molecule_id = mol2.idx), 1) == 1
        @test size(_residues(sys, molecule_id = nothing), 1) == 2
        @test size(_residues(sys, chain_id = -1), 1) == 0
        @test size(_residues(sys, chain_id = chain.idx), 1) == 1
        @test size(_residues(sys, chain_id = chain2.idx), 1) == 1
        @test size(_residues(sys, chain_id = nothing), 1) == 2
        @test size(_residues(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(_residues(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(_residues(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(_residues(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(_residues(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(_residues(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # residues_df
        df = residues_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(ResidueTuple)))
        @test copy(df[1, :]) isa ResidueTuple
        @test size(residues_df(sys), 1) == 2
        @test size(residues_df(sys, molecule_id = -1), 1) == 0
        @test size(residues_df(sys, molecule_id = mol.idx), 1) == 1
        @test size(residues_df(sys, molecule_id = mol2.idx), 1) == 1
        @test size(residues_df(sys, molecule_id = nothing), 1) == 2
        @test size(residues_df(sys, chain_id = -1), 1) == 0
        @test size(residues_df(sys, chain_id = chain.idx), 1) == 1
        @test size(residues_df(sys, chain_id = chain2.idx), 1) == 1
        @test size(residues_df(sys, chain_id = nothing), 1) == 2
        @test size(residues_df(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(residues_df(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(residues_df(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(residues_df(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(residues_df(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(residues_df(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # residues
        fv = residues(sys)
        @test fv isa Vector{Residue{T}}
        @test length(fv) == 2
        @test length(residues(sys)) == 2
        @test length(residues(sys, molecule_id = -1)) == 0
        @test length(residues(sys, molecule_id = mol.idx)) == 1
        @test length(residues(sys, molecule_id = mol2.idx)) == 1
        @test length(residues(sys, molecule_id = nothing)) == 2
        @test length(residues(sys, chain_id = -1)) == 0
        @test length(residues(sys, chain_id = chain.idx)) == 1
        @test length(residues(sys, chain_id = chain2.idx)) == 1
        @test length(residues(sys, chain_id = nothing)) == 2
        @test length(residues(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(residues(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(residues(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(residues(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(residues(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(residues(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # eachresidue
        @test length(eachresidue(sys)) == 2
        @test length(eachresidue(sys, molecule_id = -1)) == 0
        @test length(eachresidue(sys, molecule_id = mol.idx)) == 1
        @test length(eachresidue(sys, molecule_id = mol2.idx)) == 1
        @test length(eachresidue(sys, molecule_id = nothing)) == 2
        @test length(eachresidue(sys, chain_id = -1)) == 0
        @test length(eachresidue(sys, chain_id = chain.idx)) == 1
        @test length(eachresidue(sys, chain_id = chain2.idx)) == 1
        @test length(eachresidue(sys, chain_id = nothing)) == 2
        @test length(eachresidue(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(eachresidue(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(eachresidue(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(eachresidue(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(eachresidue(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(eachresidue(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # nresidues + push!
        @test nresidues(sys) isa Int
        @test nresidues(sys) == 2
        @test nresidues(sys, molecule_id = -1) == 0
        @test nresidues(sys, molecule_id = mol.idx) == 1
        @test nresidues(sys, molecule_id = mol2.idx) == 1
        @test nresidues(sys, molecule_id = nothing) == 2
        @test nresidues(sys, chain_id = -1) == 0
        @test nresidues(sys, chain_id = chain.idx) == 1
        @test nresidues(sys, chain_id = chain2.idx) == 1
        @test nresidues(sys, chain_id = nothing) == 2
        @test nresidues(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nresidues(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nresidues(sys, molecule_id = mol.idx, chain_id = chain.idx) == 1
        @test nresidues(sys, molecule_id = mol.idx, chain_id = nothing) == 1
        @test nresidues(sys, molecule_id = nothing, chain_id = chain.idx) == 1
        @test nresidues(sys, molecule_id = nothing, chain_id = nothing) == 2

        @test push!(chain, ResidueTuple(1, AminoAcid('A'))) === chain
        @test nresidues(sys) isa Int
        @test nresidues(sys) == 3
        @test nresidues(sys, molecule_id = -1) == 0
        @test nresidues(sys, molecule_id = mol.idx) == 2
        @test nresidues(sys, molecule_id = mol2.idx) == 1
        @test nresidues(sys, molecule_id = nothing) == 3
        @test nresidues(sys, chain_id = -1) == 0
        @test nresidues(sys, chain_id = chain.idx) == 2
        @test nresidues(sys, chain_id = chain2.idx) == 1
        @test nresidues(sys, chain_id = nothing) == 3
        @test nresidues(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nresidues(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nresidues(sys, molecule_id = mol.idx, chain_id = chain.idx) == 2
        @test nresidues(sys, molecule_id = mol.idx, chain_id = nothing) == 2
        @test nresidues(sys, molecule_id = nothing, chain_id = chain.idx) == 2
        @test nresidues(sys, molecule_id = nothing, chain_id = nothing) == 3

        # chain/molecule residues
        mol3 = Molecule(sys)
        @test size(_residues(mol3), 1) == 0
        @test _residues(mol3) == _residues(sys, molecule_id = mol3.idx)
        @test size(residues_df(mol3), 1) == 0
        @test residues_df(mol3) == residues_df(sys, molecule_id = mol3.idx)
        @test size(residues(mol3), 1) == 0
        @test residues(mol3) == residues(sys, molecule_id = mol3.idx)
        @test length(eachresidue(mol3)) == 0
        @test length(eachresidue(mol3)) == length(eachresidue(sys, molecule_id = mol3.idx))
        @test nresidues(mol3) == 0
        @test nresidues(mol3) == nresidues(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(_residues(chain3), 1) == 0
        @test _residues(chain3) == _residues(sys, chain_id = chain3.idx)
        @test size(residues_df(chain3), 1) == 0
        @test residues_df(chain3) == residues_df(sys, chain_id = chain3.idx)
        @test size(residues(chain3), 1) == 0
        @test residues(chain3) == residues(sys, chain_id = chain3.idx)
        @test length(eachresidue(chain3)) == 0
        @test length(eachresidue(chain3)) == length(eachresidue(sys, chain_id = chain3.idx))
        @test nresidues(chain3) == 0
        @test nresidues(chain3) == nresidues(sys, chain_id = chain3.idx)

        push!(chain3, ResidueTuple(1, AminoAcid('A')))
        @test size(_residues(mol3), 1) == 1
        @test size(_residues(mol3)) == size(_residues(sys, molecule_id = mol3.idx))
        @test size(residues_df(mol3), 1) == 1
        @test size(residues_df(mol3)) == size(residues_df(sys, molecule_id = mol3.idx))
        @test size(residues(mol3), 1) == 1
        @test residues(mol3) == residues(sys, molecule_id = mol3.idx)
        @test length(eachresidue(mol3)) == 1
        @test length(eachresidue(mol3)) == length(eachresidue(sys, molecule_id = mol3.idx))
        @test nresidues(mol3) == 1
        @test nresidues(mol3) == nresidues(sys, molecule_id = mol3.idx)

        @test size(_residues(chain3), 1) == 1
        @test size(_residues(chain3)) == size(_residues(sys, chain_id = chain3.idx))
        @test size(residues_df(chain3), 1) == 1
        @test size(residues_df(chain3)) == size(residues_df(sys, chain_id = chain3.idx))
        @test size(residues(chain3), 1) == 1
        @test residues(chain3) == residues(sys, chain_id = chain3.idx)
        @test length(eachresidue(chain3)) == 1
        @test length(eachresidue(chain3)) == length(eachresidue(sys, chain_id = chain3.idx))
        @test nresidues(chain3) == 1
        @test nresidues(chain3) == nresidues(sys, chain_id = chain3.idx)

        # residue atoms
        @test length(collect(_atoms(res))) == 0
        @test Tables.materializer(AtomTable{T})(_atoms(res)) ==
            Tables.materializer(AtomTable{T})(_atoms(sys, residue_id = res.idx))
        @test size(atoms_df(res), 1) == 0
        @test atoms_df(res) == atoms_df(sys, residue_id = res.idx)
        @test length(atoms(res)) == 0
        @test atoms(res) == atoms(sys, residue_id = res.idx)
        @test length(collect(eachatom(res))) == 0
        @test length(collect(eachatom(res))) == length(collect(eachatom(sys, residue_id = res.idx)))
        @test natoms(res) == 0
        @test natoms(res) == natoms(sys, residue_id = res.idx)

        @test push!(res, AtomTuple{T}(1, Elements.H)) === res
        @test length(collect(_atoms(res))) == 1
        @test size(atoms_df(res), 1) == 1
        @test size(atoms_df(res)) == size(atoms_df(sys, residue_id = res.idx))
        @test length(atoms(res)) == 1
        @test atoms(res) == atoms(sys, residue_id = res.idx)
        @test length(collect(eachatom(res))) == 1
        @test length(collect(eachatom(res))) == length(collect(eachatom(sys, residue_id = res.idx)))
        @test natoms(res) == 1
        @test natoms(res) == natoms(sys, residue_id = res.idx)

        # residue bonds
        @test size(_bonds(res), 1) == 0
        @test _bonds(res) == _bonds(sys, residue_id = res.idx)
        @test size(bonds_df(res), 1) == 0
        @test bonds_df(res) == bonds_df(sys, residue_id = res.idx)
        @test length(bonds(res)) == 0
        @test bonds(res) == bonds(sys, residue_id = res.idx)
        @test length(eachbond(res)) == 0
        @test length(eachbond(res)) == length(eachbond(sys, residue_id = res.idx))
        @test nbonds(res) == 0
        @test nbonds(res) == nbonds(sys, residue_id = res.idx)

        @test push!(res, BondTuple(
            Atom(res, 1, Elements.H).idx,
            Atom(res, 2, Elements.C).idx,
            BondOrder.Single
        )) === res
        @test size(_bonds(res), 1) == 1
        @test size(_bonds(res)) == size(_bonds(sys, residue_id = res.idx))
        @test size(bonds_df(res), 1) == 1
        @test size(bonds_df(res)) == size(bonds_df(sys, residue_id = res.idx))
        @test length(bonds(res)) == 1
        @test bonds(res) == bonds(sys, residue_id = res.idx)
        @test length(eachbond(res)) == 1
        @test length(eachbond(res)) == length(eachbond(sys, residue_id = res.idx))
        @test nbonds(res) == 1
        @test nbonds(res) == nbonds(sys, residue_id = res.idx)
    end
end
