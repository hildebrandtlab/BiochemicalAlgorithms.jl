
using BiochemicalAlgorithms
using CSV
using DataFramesMeta
using DelimitedFiles

# export atomtype_gaff, build_connectivity_matrix

function atomtype_gaff(num::Int64, mol::Molecule)
    element_gaff = Int((filter(:number => n -> n == num, mol.atoms).element[1]))
    
    ### Get Dataframe from csv file or from df_maker function ###
    df_ATD = CSV.read("atomtype_gff.csv", DataFrame)
    #df_ATD = df_maker("src/mappings/antechamber/ATOMTYPE_GFF.DEF")
    ### To Do: df_maker erstellt mit Typ Symbol. Schlecht für späteres Filtern mit Typ Int
    ###
    
    ### loop see below doesn't work (??) ###
    #for col in eachcol(df_ATD)
    #    replace(col, missing => -1)
       # df_ATD[Str(name)] = replace(df_ATD[Str(name)], missing => '*')
    #end
    ###
    
    #df_ATD.atomic_number = replace(df_ATD.atomic_number, missing => -1)
    #df_ATD.num_single_bonds = replace(df_ATD.num_neighbors, missing => -1)
    df_ATD = filter(:atomic_number => x -> x == element_gaff, df_ATD)
    

    bond_matrix = build_connectivity_matrix(mol)
    connected_atoms = 0
    for i = (1:lastindex(bond_matrix[num]))      
        if bond_matrix[num][i] >= 1
            connected_atoms += 1
        end
    end

    df_ATD = filter(:num_neighbors => n -> n == connected_atoms, df_ATD)

end

function build_connectivity_matrix(mol::Molecule)
    # bonds_gaff = filter(:a1 => n -> n == num, mol.bonds) #.a2
    bond_matrix = [[0 for y = (1:nrow(mol.bonds)+1)] for z = 1:nrow(mol.bonds)+1]
    # bond_matrix = Array{Int64}(undef, nrow(mol.bonds), nrow(mol.bonds)) # [[] for _ = (1:nrow(mol.bonds))]
    for i = (1:nrow(mol.bonds))
        if Int(mol.bonds.order[i]) == 1
            bond_matrix[mol.bonds.a1[i]][mol.bonds.a2[i]] = 1
            bond_matrix[mol.bonds.a2[i]][mol.bonds.a1[i]] = 1
        elseif Int(mol.bonds.order[i]) == 2
            bond_matrix[mol.bonds.a1[i]][mol.bonds.a2[i]] = 2
            bond_matrix[mol.bonds.a2[i]][mol.bonds.a1[i]] = 2
        elseif Int(mol.bonds.order[i]) == 3
            bond_matrix[mol.bonds.a1[i]][mol.bonds.a2[i]] = 3
            bond_matrix[mol.bonds.a2[i]][mol.bonds.a1[i]] = 3
        elseif Int(mol.bonds.order[i]) == 4
            bond_matrix[mol.bonds.a1[i]][mol.bonds.a2[i]] = 4
            bond_matrix[mol.bonds.a2[i]][mol.bonds.a1[i]] = 4
        end
    end
    return bond_matrix
end
