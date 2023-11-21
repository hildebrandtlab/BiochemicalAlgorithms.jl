@testitem "Protein" begin
    using BiochemicalAlgorithms: _proteins, _atoms, _bonds
    using DataFrames

    for T in [Float32, Float64]
        sys = System{T}()

        # constructors + parent
        prot = Protein(sys)
        @test prot isa Protein{T}
        @test parent(prot) === sys
        @test parent_system(prot) === sys

        if T == Float32
            prot_ds = Protein()
            @test parent(prot_ds) === default_system()
            @test parent_system(prot_ds) === default_system()

            Protein("something", Properties(:a => "b"), Flags([:A]))
        end

        prot2 = Protein(sys, "something", Properties(:a => 1), Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(prot, :_row)) == 4

        # getproperty
        @test prot.idx isa Int
        @test prot.name isa String
        @test prot.name == ""
        @test prot.properties isa Properties
        @test prot.properties == Properties()
        @test prot.flags isa Flags
        @test prot.flags == Flags()

        @test prot._sys isa System{T}
        @test prot._row isa DataFrameRow

        @test prot2.name == "something"
        @test prot2.properties == Properties(:a => 1)
        @test prot2.flags == Flags([:A, :B])

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

        # protein_by_idx
        @test_throws KeyError protein_by_idx(sys, -1)
        @test protein_by_idx(sys, prot.idx) isa Protein{T}
        @test protein_by_idx(sys, prot.idx) == prot

        # _proteins
        df = _proteins(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(ProteinTuple)))
        @test copy(df[1, :]) isa ProteinTuple

        # proteins_df
        df = proteins_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(ProteinTuple)))
        @test copy(df[1, :]) isa ProteinTuple

        # proteins
        mv = proteins(sys)
        @test mv isa Vector{Protein{T}}
        @test length(mv) == 2

        # eachprotein
        @test first(eachprotein(sys)) isa Protein{T}
        @test length(eachprotein(sys)) == 2

        # nproteins
        @test nproteins(sys) isa Int
        @test nproteins(sys) == 2

        # protein atoms
        @test size(_atoms(prot), 1) == 0
        @test _atoms(prot) == _atoms(sys, molecule_id = prot.idx)
        @test size(atoms_df(prot), 1) == 0
        @test atoms_df(prot) == atoms_df(sys, molecule_id = prot.idx)
        @test length(atoms(prot)) == 0
        @test atoms(prot) == atoms(sys, molecule_id = prot.idx)
        @test length(eachatom(prot)) == 0
        @test length(eachatom(prot)) == length(eachatom(sys, molecule_id = prot.idx))
        @test natoms(prot) == 0
        @test natoms(prot) == natoms(sys, molecule_id = prot.idx)

        @test push!(prot, AtomTuple{T}(1, Elements.H)) === prot
        @test size(_atoms(prot), 1) == 1
        @test size(_atoms(prot)) == size(_atoms(sys, molecule_id = prot.idx))
        @test size(atoms_df(prot), 1) == 1
        @test size(atoms_df(prot)) == size(atoms_df(sys, molecule_id = prot.idx))
        @test length(atoms(prot)) == 1
        @test atoms(prot) == atoms(sys, molecule_id = prot.idx)
        @test length(eachatom(prot)) == 1
        @test length(eachatom(prot)) == length(eachatom(sys, molecule_id = prot.idx))
        @test natoms(prot) == 1
        @test natoms(prot) == natoms(sys, molecule_id = prot.idx)

        for atom in eachatom(prot)
            @test parent_protein(atom) === prot
        end
        @test parent_protein(Atom(prot, 2, Elements.C)) === prot

        # protein bonds
        @test size(_bonds(prot), 1) == 0
        @test _bonds(prot) == _bonds(sys, molecule_id = prot.idx)
        @test size(bonds_df(prot), 1) == 0
        @test bonds_df(prot) == bonds_df(sys, molecule_id = prot.idx)
        @test length(bonds(prot)) == 0
        @test bonds(prot) == bonds(sys, molecule_id = prot.idx)
        @test length(eachbond(prot)) == 0
        @test length(eachbond(prot)) == length(eachbond(sys, molecule_id = prot.idx))
        @test nbonds(prot) == 0
        @test nbonds(prot) == nbonds(sys, molecule_id = prot.idx)

        @test push!(prot, BondTuple(
            Atom(prot, 1, Elements.H).idx,
            Atom(prot, 2, Elements.C).idx,
            BondOrder.Single)
        ) === prot
        @test size(_bonds(prot), 1) == 1
        @test size(_bonds(prot)) == size(_bonds(sys, molecule_id = prot.idx))
        @test size(bonds_df(prot), 1) == 1
        @test size(bonds_df(prot)) == size(bonds_df(sys, molecule_id = prot.idx))
        @test length(bonds(prot)) == 1
        @test bonds(prot) == bonds(sys, molecule_id = prot.idx)
        @test length(eachbond(prot)) == 1
        @test length(eachbond(prot)) == length(eachbond(sys, molecule_id = prot.idx))
        @test nbonds(prot) == 1
        @test nbonds(prot) == nbonds(sys, molecule_id = prot.idx)
    end
end
