
using BiochemicalAlgorithms: Molecule, Bond, Atom
using DataFrames, DelimitedFiles, CSV

export Atomtype, make_dict_from_dat


function atomtype_gaff(num::UInt128, mol::Molecule)
    atom_gaff = mol.atoms[num]
    bonds_gaff = subset(mol.bonds, aid1 == num)
    
    for i = 1:length(max(Molecule.bonds.a1))
    end
end
