export
    Protein,
    ProteinTable,
    protein_by_idx,
    proteins,
    nproteins,
    parent_protein

# TODO implement rules to distinguish from molecules
const ProteinTable{T} = MoleculeTable{T}
const Protein{T} = Molecule{T}
protein_by_idx(sys::System, idx::Int) = molecule_by_idx(sys, idx)
proteins(sys::System) = molecules(sys)
nproteins(sys::System) = nmolecules(sys)
parent_protein(ac::Union{Atom, Chain, Fragment, Nucleotide, Residue}) = parent_molecule(ac)
