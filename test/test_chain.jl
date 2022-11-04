using DataFrames

@testset "PDBChain" begin

    # test default constructor
    c1 = PDBChain()
    @test c1 isa PDBChain
    @test c1.name == ""
    @test c1.fragments isa DataFrame
    @test fragments(c1) == c1.fragments

    # test costume constructor
    f = DataFrame(number=[1], name=["fra"], chain_id=["A"])
    c2 = PDBChain("Chain A", f)
    @test c2.name == "Chain A"
    @test fragments(c2) isa DataFrame

    names = keys(eachcol(fragments(c2)))
    types = eltype.(eachcol(fragments(c2)))
    for i in eachindex(names)
        @test names[i] == keys(eachcol(f))[i]
        @test types[i] == eltype.(eachcol(f))[i] 
    end
end

@testset "ProteinChain" begin

    # test default constructor
    c3 = ProteinChain()
    @test c3 isa ProteinChain
    @test c3.name == ""
    @test c3.residues isa DataFrame
    @test fragments(c3) == c3.residues

    # test costume constructor
    f = DataFrame(number=[1], name=["fra"], chain_id=["A"])
    c4 = ProteinChain("Chain A", f)
    @test c4.name == "Chain A"
    @test fragments(c4) isa DataFrame
   
    names =  keys(eachcol(fragments(c4)))
    types = eltype.(eachcol(fragments(c4)))
    
    for i in eachindex(names)
        @test names[i] == keys(eachcol(f))[i]
        @test types[i] == eltype.(eachcol(f))[i] 
    end
end

@testset "HeteroAtomChain" begin
    # test default constructor
    c5 = HeteroAtomChain()
    @test c5 isa HeteroAtomChain
    @test c5.name == ""
    @test c5.fragments isa DataFrame
    @test fragments(c5) == c5.fragments

    # test costume constructor
    f = DataFrame(number=[1], name=["fra"], chain_id=["A"])
    c6 = HeteroAtomChain("Chain A", f)
    @test c6.name == "Chain A"
    @test fragments(c6) isa DataFrame
   
    names =  keys(eachcol(fragments(c6)))
    types = eltype.(eachcol(fragments(c6)))
    
    for i in eachindex(names)
        @test names[i] == keys(eachcol(f))[i]
        @test types[i] == eltype.(eachcol(f))[i] 
    end
end

@testset "NucleotideChain" begin
    # test default constructor
    c7 = NucleotideChain()
    @test c7 isa NucleotideChain
    @test c7.name == ""
    @test c7.nucleotides isa DataFrame
    @test fragments(c7) == c7.nucleotides

    # test costume constructor
    f = DataFrame(number=[1], name=["fra"], chain_id=["A"])
    c8 = NucleotideChain("Chain A", f)
    @test c8.name == "Chain A"
    @test fragments(c8) isa DataFrame
   
    names =  keys(eachcol(fragments(c8)))
    types = eltype.(eachcol(fragments(c8)))
    
    for i in eachindex(names)
        @test names[i] == keys(eachcol(f))[i]
        @test types[i] == eltype.(eachcol(f))[i] 
    end
end
