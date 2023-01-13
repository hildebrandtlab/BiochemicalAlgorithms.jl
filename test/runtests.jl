using BiochemicalAlgorithms
using Test

@testset "BiochemicalAlgorithms.jl" begin
   
    @testset "Fileformats" begin include("test_fileformats.jl") end
    @testset "Core" begin 
        include("test_amino_acid.jl")
        include("test_atom.jl")
        include("test_bond.jl")
        include("test_chain.jl")
        include("test_element.jl")
        include("test_fragment.jl")
        include("test_molecule.jl")
        include("test_nucleotide.jl")
        include("test_pdb_atom.jl")
        include("test_pdb_molecule.jl")
        include("test_protein.jl")
        include("test_residue.jl")
        include("test_types.jl")
    end
end
