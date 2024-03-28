export
    FragmentVariantType,
    FragmentVariant,
    MoleculeVariantType,
    MoleculeVariant

@enumx MoleculeVariant begin
    None = 1
    Protein = 2
end

const MoleculeVariantType = MoleculeVariant.T

@enumx FragmentVariant begin
    None = 1
    Residue = 2
    Nucleotide = 3
end

const FragmentVariantType = FragmentVariant.T
