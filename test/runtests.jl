using BiochemicalAlgorithms
using Test
using DataFrames

@testset verbose=true "BiochemicalAlgorithms.jl" begin
   
    @testset "Core" begin 
        include("core/test_amino_acid.jl")
        include("core/test_bond_order.jl")
        include("core/test_atom.jl")
        include("core/test_bond.jl")
        include("core/test_element.jl")
        include("core/test_fragment.jl")
        include("core/test_molecule.jl")
        include("core/test_nucleotide.jl")
        include("core/test_protein.jl")
        include("core/test_residue.jl")
        include("core/test_types.jl")
        include("core/test_tuples.jl")
    end

    @testset "Fileformats" begin 
        include("fileformats/test_pdb.jl")
        include("fileformats/test_pubchem_json.jl")
        include("fileformats/test_sdfile.jl")
    end
    
    @testset "Preprocessing" begin
        include("preprocessing/test_fragmentdb.jl")
    end

    @testset "Substructures" begin
        include("test_substructure.jl")
    end
end
