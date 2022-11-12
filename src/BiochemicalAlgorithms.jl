module BiochemicalAlgorithms

include("core/types.jl")
include("core/element.jl")
include("core/amino_acid.jl")
include("core/atom.jl")
include("core/PDB_atom.jl")
include("core/bond.jl")
include("core/chain.jl")
include("core/molecule.jl")
include("core/pdb_molecule.jl")
include("core/residue.jl")
include("core/fragment.jl")
include("core/nucleotide.jl")
include("core/protein.jl")
include("core/atomtyping.jl")
include("core/atomtyping_prep.jl")
include("core/def_bond.jl")

# testing script to load all molecules from folder
include("validation/atomtypes_dict_script.jl")
include("validation/validation_script.jl")

module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/PDB.jl")
include("fileformats/mol2.jl")

include("mappings/atom_bijection.jl")
include("mappings/rigid_mapping.jl")

using .PubChem

export load_pubchem_json

end
