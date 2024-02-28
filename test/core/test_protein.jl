@testitem "Protein" begin
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
        @test prot._row isa BiochemicalAlgorithms._MoleculeTableRow

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

        # proteins
        mv = proteins(sys)
        @test mv isa ProteinTable{T}
        @test length(mv) == 2

        # protein atoms
        @test length(atoms(prot)) == 0
        @test atoms(prot) == atoms(sys, molecule_id = prot.idx)
        @test natoms(prot) == 0
        @test natoms(prot) == natoms(sys, molecule_id = prot.idx)

        @test push!(prot, AtomTuple{T}(1, Elements.H)) === prot
        @test length(atoms(prot)) == 1
        @test atoms(prot) == atoms(sys, molecule_id = prot.idx)
        @test natoms(prot) == 1
        @test natoms(prot) == natoms(sys, molecule_id = prot.idx)

        for atom in atoms(prot)
            @test parent_protein(atom) === prot
        end
        @test parent_protein(Atom(prot, 2, Elements.C)) === prot

        # protein bonds
        @test length(bonds(prot)) == 0
        @test bonds(prot) == bonds(sys, molecule_id = prot.idx)
        @test nbonds(prot) == 0
        @test nbonds(prot) == nbonds(sys, molecule_id = prot.idx)

        @test push!(prot, BondTuple(
            Atom(prot, 1, Elements.H).idx,
            Atom(prot, 2, Elements.C).idx,
            BondOrder.Single)
        ) === prot
        @test length(bonds(prot)) == 1
        @test bonds(prot) == bonds(sys, molecule_id = prot.idx)
        @test nbonds(prot) == 1
        @test nbonds(prot) == nbonds(sys, molecule_id = prot.idx)
    end
end
