module BiochemicalAlgorithms

using AutoHashEquals
using CellListMap
using DataFrames
using DataStructures
using DocStringExtensions
using Graphs, GraphDataFrameBridge
using PrettyTables
using StaticArrays
using Tables, TableOperations
using Unitful, UnitfulAtomic

include("core/exceptions.jl")
include("core/constants.jl")
include("core/types.jl")
include("core/element.jl")
include("core/amino_acid.jl")
include("core/bond_order.jl")
include("core/system_component.jl")
include("core/atom_container.jl")
include("core/system_atoms.jl")
include("core/system_bonds.jl")
include("core/system_molecules.jl")
include("core/system_chains.jl")
include("core/system_fragments.jl")
include("core/system_nucleotides.jl")
include("core/system_residues.jl")
include("core/system.jl")
include("core/system_tables.jl")
include("core/bond.jl")
include("core/atom.jl")
include("core/molecule.jl")
include("core/chain.jl")
include("core/fragment.jl")
include("core/nucleotide.jl")
include("core/residue.jl")
include("core/moleculargraph_wrapper.jl")

include("substructures/substructure.jl")
include("substructures/smarts.jl")
include("substructures/sssr.jl")

include("fileformats/ball_ini_file.jl")
module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/pdb.jl")
include("fileformats/sdfile.jl")

include("mappings/atom_bijection.jl")
include("mappings/rigid_mapping.jl")

include("forcefields/common/forcefield_parameters.jl")
include("forcefields/common/atomtype_template.jl")
include("forcefields/common/forcefield.jl")
include("forcefields/common/stretch_component.jl")
include("forcefields/common/bend_component.jl")
include("forcefields/common/torsion_component.jl")
include("forcefields/common/nonbonded_component.jl")

include("forcefields/AMBER/amberff_parameters.jl")
include("forcefields/AMBER/amberff.jl")

include("forcefields/MMFF94/mmff94_parameters.jl")

include("preprocessing/fragmentDB.jl")
include("preprocessing/normalize_names.jl")
include("preprocessing/build_bonds.jl")
include("preprocessing/add_hydrogens.jl")
include("preprocessing/reconstruct_fragments.jl")

using .PubChem

export load_pubchem_json, ball_data_path

ball_data_path(parts...) = normpath(joinpath(@__DIR__, "..", "data", parts...))

end
