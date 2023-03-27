export Protein, protein_by_idx, proteins, proteins_df, eachprotein, nproteins

# TODO implement rules to distinguish from molecules
const Protein{T} = Molecule{T}
_proteins(sys::System) = _molecules(sys)
protein_by_idx(sys::System, idx::Int) = molecule_by_idx(sys, idx)
proteins(sys::System) = molecules(sys)
proteins_df(sys::System) = molecules_df(sys)
eachprotein(sys::System) = eachmolecule(sys)
nproteins(sys::System) = nmolecules(sys)
