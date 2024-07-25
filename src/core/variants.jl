export
    FragmentVariantType,
    FragmentVariant,
    MoleculeVariantType,
    MoleculeVariant

"""
    @enumx MoleculeVariant begin
        None = 1
        Protein = 2
    end

Enum representing variants of molecules

# Example
```jldoctest
julia> prot = Protein(System())
Molecule{Float32}: (idx = 1, name = "")

julia> isprotein(prot)
true

julia> prot.variant == MoleculeVariant.Protein
true
```
"""
@enumx MoleculeVariant begin
    None = 1
    Protein = 2
end

"""
    const MoleculeVariantType = MoleculeVariant.T

Type of `MoleculeVariant` enumerators
"""
const MoleculeVariantType = MoleculeVariant.T

"""
    @enumx FragmentVariant begin
        None = 1
        Residue = 2
        Nucleotide = 3
    end

Enum representing variants of fragments

# Example
```jldoctest
julia> res = Residue(Chain(Molecule(System())), 1; name = "ALA")
Fragment{Float32}: (idx = 3, number = 1, name = "ALA")

julia> isresidue(res)
true

julia> res.variant == FragmentVariant.Residue
true
```
"""
@enumx FragmentVariant begin
    None = 1
    Residue = 2
    Nucleotide = 3
end

"""
    const FragmentVariantType = FragmentVariant.T

Type of `FragmentVariant` enumerators
"""
const FragmentVariantType = FragmentVariant.T
