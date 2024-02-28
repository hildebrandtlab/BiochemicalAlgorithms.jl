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
        @test length(getfield(res, :_row)) == 5

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
        @test res._row isa BiochemicalAlgorithms._ResidueTableRow
        
        @test res.molecule_id isa Int
        @test res.molecule_id == mol.idx
        @test res.chain_id isa Int
        @test res.chain_id == chain.idx

        @test res2.number == 1
        @test res2.type == AminoAcid('D')
        @test res2.properties == Properties(:a => 1)
        @test res2.flags == Flags([:A, :B])
        @test res2.molecule_id == mol2.idx
        @test res2.chain_id == chain2.idx

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
        res3.molecule_id = 999
        @test res3.molecule_id == 999
        res3.chain_id = 998
        @test res3.chain_id == 998

        # residue_by_idx
        @test_throws KeyError residue_by_idx(sys, -1)
        @test residue_by_idx(sys, res.idx) isa Residue{T}
        @test residue_by_idx(sys, res.idx) == res

        # residues
        fv = residues(sys)
        @test fv isa ResidueTable{T}
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
        @test size(residues(mol3), 1) == 0
        @test residues(mol3) == residues(sys, molecule_id = mol3.idx)
        @test nresidues(mol3) == 0
        @test nresidues(mol3) == nresidues(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(residues(chain3), 1) == 0
        @test residues(chain3) == residues(sys, chain_id = chain3.idx)
        @test nresidues(chain3) == 0
        @test nresidues(chain3) == nresidues(sys, chain_id = chain3.idx)

        push!(chain3, ResidueTuple(1, AminoAcid('A')))
        @test size(residues(mol3), 1) == 1
        @test residues(mol3) == residues(sys, molecule_id = mol3.idx)
        @test nresidues(mol3) == 1
        @test nresidues(mol3) == nresidues(sys, molecule_id = mol3.idx)

        @test size(residues(chain3), 1) == 1
        @test residues(chain3) == residues(sys, chain_id = chain3.idx)
        @test nresidues(chain3) == 1
        @test nresidues(chain3) == nresidues(sys, chain_id = chain3.idx)

        # residue atoms
        @test length(atoms(res)) == 0
        @test atoms(res) == atoms(sys, residue_id = res.idx)
        @test natoms(res) == 0
        @test natoms(res) == natoms(sys, residue_id = res.idx)

        @test push!(res, AtomTuple{T}(1, Elements.H)) === res
        @test length(atoms(res)) == 1
        @test atoms(res) == atoms(sys, residue_id = res.idx)
        @test natoms(res) == 1
        @test natoms(res) == natoms(sys, residue_id = res.idx)

        # residue bonds
        @test length(bonds(res)) == 0
        @test bonds(res) == bonds(sys, residue_id = res.idx)
        @test nbonds(res) == 0
        @test nbonds(res) == nbonds(sys, residue_id = res.idx)

        @test push!(res, BondTuple(
            Atom(res, 1, Elements.H).idx,
            Atom(res, 2, Elements.C).idx,
            BondOrder.Single
        )) === res
        @test length(bonds(res)) == 1
        @test bonds(res) == bonds(sys, residue_id = res.idx)
        @test nbonds(res) == 1
        @test nbonds(res) == nbonds(sys, residue_id = res.idx)
    end
end
