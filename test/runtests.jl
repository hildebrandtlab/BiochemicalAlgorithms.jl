using BiochemicalAlgorithms
using Test

@testset "BiochemicalAlgorithms.jl" begin
   
    @testset "Fileformats" begin include("test_fileformats.jl") end
    @testset "Core" begin 
            include("test_amino_acid.jl") 
        end
end
