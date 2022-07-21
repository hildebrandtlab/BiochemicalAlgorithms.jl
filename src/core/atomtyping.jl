
using BiochemicalAlgorithms: Molecule, Bond, Atom
using DataFrames

export Atomtype, make_dict_from_dat

function atomtype_gaff(num::UInt128)
    atom_gaff = BiochemicalAlgorithms.Molecule.atoms[num]
    # bonds_gaff = subset(groupby(BiochemicalAlgorithms.Molecule.bonds, :order), :aid2 => aid1 -> aid1 = num)
    bonds_gaff = subset(BiochemicalAlgorithms.Molecule.bonds, aid1 == num)
    
    for i = 1:length(max(Molecule.bonds.a1))
    
end

function make_dict_from_dat(mapfile::AbstractString, firstline::UInt16, lastline::UInt16)
    dat_file = open(mapfile)
    Dat_dict = Dict{String, UInt16}()
    for (i, line) in enumerate(eachline(dat_file))
        if (i >= firstline & i <= lastline)
            Dat_dict(line[1:2] => i)
        end
    end
end