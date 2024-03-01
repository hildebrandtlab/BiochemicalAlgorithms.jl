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
        @test length(getfield(nuc, :_row)) == 5

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

        @test nuc._sys isa System{T}
        @test nuc._row isa BiochemicalAlgorithms._NucleotideTableRow
        
        @test nuc.molecule_id isa Int
        @test nuc.molecule_id == mol.idx
        @test nuc.chain_id isa Int
        @test nuc.chain_id == chain.idx

        @test nuc2.number == 1
        @test nuc2.name == "something"
        @test nuc2.properties == Properties(:a => 1)
        @test nuc2.flags == Flags([:A, :B])
        @test nuc2.molecule_id == mol2.idx
        @test nuc2.chain_id == chain2.idx

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
        nuc3.molecule_id = 999
        @test nuc3.molecule_id == 999
        nuc3.chain_id = 998
        @test nuc3.chain_id == 998

        # nucleotide_by_idx
        @test_throws KeyError nucleotide_by_idx(sys, -1)
        @test nucleotide_by_idx(sys, nuc.idx) isa Nucleotide{T}
        @test nucleotide_by_idx(sys, nuc.idx) == nuc

        # nucleotides
        fv = nucleotides(sys)
        @test fv isa NucleotideTable{T}
        @test length(fv) == 2
        @test length(nucleotides(sys)) == 2
        @test length(nucleotides(sys, molecule_id = -1)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx)) == 1
        @test length(nucleotides(sys, molecule_id = mol2.idx)) == 1
        @test length(nucleotides(sys, molecule_id = nothing)) == 2
        @test length(nucleotides(sys, chain_id = -1)) == 0
        @test length(nucleotides(sys, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, chain_id = chain2.idx)) == 1
        @test length(nucleotides(sys, chain_id = nothing)) == 2
        @test length(nucleotides(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(nucleotides(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # nnucleotides
        @test nnucleotides(sys) isa Int
        @test nnucleotides(sys) == 2
        @test nnucleotides(sys, molecule_id = -1) == 0
        @test nnucleotides(sys, molecule_id = mol.idx) == 1
        @test nnucleotides(sys, molecule_id = mol2.idx) == 1
        @test nnucleotides(sys, molecule_id = nothing) == 2
        @test nnucleotides(sys, chain_id = -1) == 0
        @test nnucleotides(sys, chain_id = chain.idx) == 1
        @test nnucleotides(sys, chain_id = chain2.idx) == 1
        @test nnucleotides(sys, chain_id = nothing) == 2
        @test nnucleotides(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = chain.idx) == 1
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = nothing) == 1
        @test nnucleotides(sys, molecule_id = nothing, chain_id = chain.idx) == 1
        @test nnucleotides(sys, molecule_id = nothing, chain_id = nothing) == 2

        # chain/molecule nucleotides
        mol3 = Molecule(sys)
        @test size(nucleotides(mol3), 1) == 0
        @test nucleotides(mol3) == nucleotides(sys, molecule_id = mol3.idx)
        @test nnucleotides(mol3) == 0
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(nucleotides(chain3), 1) == 0
        @test nucleotides(chain3) == nucleotides(sys, chain_id = chain3.idx)
        @test nnucleotides(chain3) == 0
        @test nnucleotides(chain3) == nnucleotides(sys, chain_id = chain3.idx)

        Nucleotide(chain3, 1)
        @test size(nucleotides(mol3), 1) == 1
        @test nucleotides(mol3) == nucleotides(sys, molecule_id = mol3.idx)
        @test nnucleotides(mol3) == 1
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_id = mol3.idx)

        @test size(nucleotides(chain3), 1) == 1
        @test nucleotides(chain3) == nucleotides(sys, chain_id = chain3.idx)
        @test nnucleotides(chain3) == 1
        @test nnucleotides(chain3) == nnucleotides(sys, chain_id = chain3.idx)

        # nucleotide atoms
        @test length(atoms(nuc)) == 0
        @test atoms(nuc) == atoms(sys, nucleotide_id = nuc.idx)
        @test natoms(nuc) == 0
        @test natoms(nuc) == natoms(sys, nucleotide_id = nuc.idx)

        @test Atom(nuc, 1, Elements.H).nucleotide_id == nuc.idx
        @test length(atoms(nuc)) == 1
        @test atoms(nuc) == atoms(sys, nucleotide_id = nuc.idx)
        @test natoms(nuc) == 1
        @test natoms(nuc) == natoms(sys, nucleotide_id = nuc.idx)

        # nucleotide bonds
        @test length(bonds(nuc)) == 0
        @test bonds(nuc) == bonds(sys, nucleotide_id = nuc.idx)
        @test nbonds(nuc) == 0
        @test nbonds(nuc) == nbonds(sys, nucleotide_id = nuc.idx)

        Bond(nuc, Atom(nuc, 1, Elements.H).idx, Atom(nuc, 2, Elements.C).idx, BondOrder.Single)
        @test length(bonds(nuc)) == 1
        @test bonds(nuc) == bonds(sys, nucleotide_id = nuc.idx)
        @test nbonds(nuc) == 1
        @test nbonds(nuc) == nbonds(sys, nucleotide_id = nuc.idx)
    end
end
