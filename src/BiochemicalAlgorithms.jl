module BiochemicalAlgorithms

using AutoHashEquals
using CellListMap
using DataFrames
using DataStructures
using DocStringExtensions
using Graphs, GraphDataFrameBridge
using Unitful, UnitfulAtomic

import PrettyTables
import StaticArrays
import Tables, TableOperations

# core definitions
include("core/exceptions.jl")
include("core/constants.jl")
include("core/types.jl")
include("core/tables.jl")
include("core/element.jl")
include("core/amino_acid.jl")
include("core/bond_order.jl")

# system
include("core/system_internals/_atom_table.jl")
include("core/system_internals/_bond_table.jl")
include("core/system_internals/_molecule_table.jl")
include("core/system_internals/_chain_table.jl")
include("core/system_internals/_fragment_table.jl")
include("core/system_internals/_nucleotide_table.jl")
include("core/system_internals/_residue_table.jl")
include("core/system.jl")

# system components
include("core/system_component_table.jl")
include("core/bond.jl")
include("core/atom.jl")
include("core/molecule.jl")
include("core/chain.jl")
include("core/fragment.jl")
include("core/nucleotide.jl")
include("core/residue.jl")

# molgraph
include("core/moleculargraph_wrapper.jl")

# substructures
include("substructures/substructure.jl")
include("substructures/smarts.jl")
include("substructures/sssr.jl")

# file formats
include("fileformats/ball_ini_file.jl")
module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/pdb.jl")
include("fileformats/sdfile.jl")

# mappings
include("mappings/atom_bijection.jl")
include("mappings/rigid_mapping.jl")

# force fields
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

# preprocessing
include("preprocessing/fragmentDB.jl")
include("preprocessing/normalize_names.jl")
include("preprocessing/build_bonds.jl")
include("preprocessing/add_hydrogens.jl")
include("preprocessing/reconstruct_fragments.jl")

using .PubChem

export load_pubchem_json, ball_data_path

ball_data_path(parts...) = normpath(joinpath(@__DIR__, "..", "data", parts...))

end
