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
include("core/system.jl")

module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/PDB.jl")

using .PubChem

export load_pubchem_json

end
