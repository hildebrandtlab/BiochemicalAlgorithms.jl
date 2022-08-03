
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs

export start_menu_atomtype, atomtype_gaff, build_graph

# export atomtype_gaff

function start_menu_atomtype()
    println("Which atomtype do you want? \n Press: 1 for GAFF,  2 for Amber")
    def_number = readline()
    def_number = parse(Int8, chomp(def_number))
    if def_number == 1
        df_ATD = CSV.read("atomtype_gff.csv", DataFrame)
        return df_ATD
    elseif def_number == 2
        df_ATD = CSV.read("atomtype_amber.csv", DataFrame)
        return df_ATD
    end
    return "Atomtyping canceled"
end


function atomtype_gaff(df_ATD::DataFrame, mol::Molecule)
    ATD_array = Array{String, 1}(undef,nrow(mol.atoms))
    for i = (1:length(ATD_array)) 
        ATD_array[i] = String(Symbol(mol.atoms.element[i]))
    end
    mol_graph = build_graph(mol)
    adj_matrix = adjacency_matrix(mol_graph)
    for num = (1:nrow(mol.atoms))
        element_gaff = Int8((filter(:number => n -> n == num, mol.atoms).element[1]))
        df_ATD_temp = filter(:atomic_number => x -> x == element_gaff, df_ATD)
        
        connected_atoms = 0
        for i = (1:nrow(mol.atoms))      
            if adj_matrix[num,i] >= 1
                connected_atoms += 1
            end
        end

        df_ATD_temp = filter(:num_neighbors => n -> n == connected_atoms, df_ATD_temp)
        
        df_ATD_temp = CES_checker(num, mol, mol_graph, df_ATD_temp, ATD_array)
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_array[num] = df_ATD_temp.type_name[1]
        end    
    end
    return ATD_array
end


function CES_checker(num::Integer, mol::Molecule, graph::SimpleWeightedGraph, df::DataFrame, ATD_array::AbstractArray)
    if !has_self_loops(graph)
        df = filter(:atomic_property => n -> !occursin("RG", n), df)
        df = filter(:atomic_property => n -> !occursin("AR", n), df)
    end
    ### To Do: check for neighbors with Graphs. neighbors(graph, num) in 
    return df
end


function build_graph(mol::Molecule)
    mol_graph = SimpleWeightedGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        if Int8(mol.bonds.order[i]) == 1
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 1)
        elseif Int8(mol.bonds.order[i]) == 2
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 2)
        elseif Int8(mol.bonds.order[i]) == 3
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 3)
        elseif Int8(mol.bonds.order[i]) == 4
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 4)
        end
    end
    return mol_graph
end

