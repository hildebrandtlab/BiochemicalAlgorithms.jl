@testset "Atom" begin
    for T in [Float32, Float64]
        at = (
            idx = 0,
            number = 1,
            element = Elements.H,
            name = "my fancy atom",
            atomtype = "heavy",
            r = Vector3{T}(1, 2, 4),
            v = Vector3{T}(1, 1, 1),
            F = Vector3{T}(0, 0, 0),
            has_velocity = true,
            has_force = false,
            properties = Dict{String, Any}()
        )::AtomTuple{T}
        sys = System{T}()

        # constructors + parent
        atom = Atom(sys, at)
        @test atom isa Atom{T}
        @test parent(atom) === sys
        T == Float32 && @test parent(Atom(at)) === default_system()

        atom2 = Atom(sys, at, frame_id = 10, molecule_id = 11, chain_id = 12, fragment_id = 13, 
            nucleotide_id = 14, residue_id = 15)

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(atom, :row)) == 17

        # getproperty (public)
        @test atom.idx isa Int
        @test atom.number isa Int
        @test atom.number == at.number
        @test atom.element isa ElementType
        @test atom.element == at.element
        @test atom.name isa String
        @test atom.name == at.name
        @test atom.atomtype isa String
        @test atom.atomtype == at.atomtype
        @test atom.r isa Vector3{T}
        @test atom.r == at.r
        @test atom.v isa Vector3{T}
        @test atom.v == at.v
        @test atom.F isa Vector3{T}
        @test atom.F == at.F
        @test atom.has_velocity isa Bool
        @test atom.has_velocity == at.has_velocity
        @test atom.has_force isa Bool
        @test atom.has_force == at.has_force
        @test atom.properties isa Properties
        @test atom.properties == at.properties

        @test atom.sys isa System{T}
        @test atom.row isa DataFrameRow

        @test_throws ErrorException atom.frame_id
        @test_throws ErrorException atom.molecule_id
        @test_throws ErrorException atom.chain_id
        @test_throws ErrorException atom.fragment_id
        @test_throws ErrorException atom.nucleotide_id
        @test_throws ErrorException atom.residue_id

        @test atom.row.frame_id isa Int
        @test atom.row.frame_id == 1
        @test ismissing(atom.row.molecule_id)
        @test ismissing(atom.row.chain_id)
        @test ismissing(atom.row.fragment_id)
        @test ismissing(atom.row.nucleotide_id)
        @test ismissing(atom.row.residue_id)

        @test atom2.row.frame_id isa Int
        @test atom2.row.frame_id == 10
        @test atom2.row.molecule_id isa Int
        @test atom2.row.molecule_id == 11
        @test atom2.row.chain_id isa Int
        @test atom2.row.chain_id == 12
        @test atom2.row.fragment_id isa Int
        @test atom2.row.fragment_id == 13
        @test atom2.row.nucleotide_id isa Int
        @test atom2.row.nucleotide_id == 14
        @test atom2.row.residue_id isa Int
        @test atom2.row.residue_id == 15

        # setproperty!
        atom.number = 42
        @test atom.number == 42
        atom.element = Elements.C
        @test atom.element == Elements.C
        atom.name = "another name"
        @test atom.name == "another name"
        atom.atomtype = "none"
        @test atom.atomtype == "none"
        atom.r = Vector3{T}(10, 20, 30)
        @test atom.r == Vector3{T}(10, 20, 30)
        atom.v = Vector3{T}(100, 200, 300)
        @test atom.v == Vector3{T}(100, 200, 300)
        atom.F = Vector3{T}(1000, 2000, 3000)
        @test atom.F == Vector3{T}(1000, 2000, 3000)
        atom.has_velocity = false
        @test !atom.has_velocity
        atom.has_force = true
        @test atom.has_force
        atom.properties = Properties("first" => "v1", "second" => 99)
        @test atom.properties["first"] == "v1"
        @test atom.properties["second"] == 99

        @test_throws ErrorException atom.frame_id = 0
        @test_throws ErrorException atom.molecule_id = 0
        @test_throws ErrorException atom.chain_id = 0
        @test_throws ErrorException atom.fragment_id = 0
        @test_throws ErrorException atom.nucleotide_id = 0
        @test_throws ErrorException atom.residue_id = 0

        # atoms_df
        df = atoms_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (1, length(fieldnames(AtomTuple{T})))
        @test copy(df[1, :]) isa AtomTuple{T}
        @test size(atoms(sys, frame_id = 1), 1) == 1
        @test size(atoms(sys, frame_id = 2), 1) == 0
        @test size(atoms(sys, frame_id = 10), 1) == 1
        @test size(atoms(sys, frame_id = nothing), 1) == 2
        @test size(atoms(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15), 1) == 1

        # eachatom
        @test length(eachatom(sys)) == 1
        @test first(eachatom(sys)) isa Atom{T}
        @test length(eachatom(sys, frame_id = 1)) == 1
        @test length(eachatom(sys, frame_id = 2)) == 0
        @test length(eachatom(sys, frame_id = 10)) == 1
        @test length(eachatom(sys, frame_id = nothing)) == 2
        @test length(eachatom(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15)) == 1

        @test length(atoms(sys)) == 1
        @test first(atoms(sys)) isa Atom{T}
        @test length(atoms(sys, frame_id = 1)) == 1
        @test length(atoms(sys, frame_id = 2)) == 0
        @test length(atoms(sys, frame_id = 10)) == 1
        @test length(atoms(sys, frame_id = nothing)) == 2
        @test length(atoms(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15)) == 1

        # natoms + push!
        @test natoms(sys) isa Int
        @test natoms(sys) == 1
        @test natoms(sys, frame_id = 1) == 1
        @test natoms(sys, frame_id = 2) == 0
        @test natoms(sys, frame_id = 10) == 1
        @test natoms(sys, frame_id = nothing) == 2
        @test natoms(sys, frame_id = nothing, molecule_id =11, chain_id = 12, fragment_id = 13,
            nucleotide_id = 14, residue_id = 15) == 1

        @test push!(sys, at) === sys
        @test natoms(sys) == 2
        @test push!(sys, at, frame_id = 100, molecule_id = 101, chain_id = 102, fragment_id = 103, 
            nucleotide_id = 104, residue_id = 105) === sys
        @test natoms(sys) == 2
        @test natoms(sys, frame_id = 100) == 1
    end
end
