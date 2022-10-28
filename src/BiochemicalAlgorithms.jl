module BiochemicalAlgorithms

include("core/types.jl")
include("core/element.jl")
include("core/amino_acid.jl")
include("core/atom.jl")
include("core/pdb_atom.jl")
include("core/bond.jl")
include("core/chain.jl")
include("core/molecule.jl")
include("core/pdb_molecule.jl")
include("core/residue.jl")
include("core/fragment.jl")
include("core/nucleotide.jl")
include("core/protein.jl")
include("core/def_bond.jl")


module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/pdb.jl")
include("fileformats/mol2.jl")

include("mappings/atom_bijection.jl")
include("mappings/rigid_mapping.jl")

using .PubChem

export load_pubchem_json

end
