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
        end
end
