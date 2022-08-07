
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs

export select_atomtyping, get_atomtype, build_graph


function select_atomtyping()
    println("Which atomtyping do you want? \n Press: 1 for GAFF,  2 for Amber")
    def_number = readline()
    def_number = parse(Int8, chomp(def_number))
    if def_number == 1
        df_ATD = CSV.read("atomtype_gff.csv", DataFrame)
        return df_ATD
    elseif def_number == 2
        df_ATD = CSV.read("atomtype_amber.csv", DataFrame)
        return df_ATD
    end
    return "No Atomtype definition chosen. cancelled!"
end


function get_atomtype(mol::Molecule, df_ATD::DataFrame)
    ATD_array = Array{String, 1}(undef,nrow(mol.atoms))
    for i = (1:length(ATD_array)) 
        ATD_array[i] = String(Symbol(mol.atoms.element[i]))
    end
    mol_graph = build_graph(mol)
    adj_matrix = adjacency_matrix(mol_graph)
    
    # Cycle detection and list
    cycle_bool = (ne(mol_graph) >= nv(mol_graph))
    chem_cycle_list = cycle_checker(mol_graph)

    for num = (1:nrow(mol.atoms))
        # filter for element in df
        element_gaff = Int8((filter(:number => n -> n == num, mol.atoms).element[1]))
        df_ATD_temp = filter(:atomic_number => x -> x == element_gaff, df_ATD)
        df_neighbors = filter(:a1 => n -> n == num, mol.bonds)
        
        #filter for connected_atoms and connected_H_atoms
        connected_atoms = 0
        connected_H_atoms = 0
        for i = (1:nrow(mol.atoms))      
            if adj_matrix[num,i] >= 1
                connected_atoms += 1
                if Int8((filter(:number => n -> n == i, mol.atoms).element[1])) == 1
                    connected_H_atoms += 1
                end
            end
        end
        
        df_ATD_temp = filter(:num_neighbors => n -> n == connected_atoms, df_ATD_temp)
        if connected_H_atoms >= 1 && (1 in df_ATD_temp.num_H_bonds || 2 in df_ATD_temp.num_H_bonds || 3 in df_ATD_temp.num_H_bonds)
            df_ATD_temp = filter(:num_H_bonds => n -> n == connected_H_atoms, df_ATD_temp)
        end

        # filter out obvious loop properties if no cycle detected
        if !cycle_bool 
            df_ATD_temp = filter(:atomic_property => n -> !occursin("RG", n), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> !occursin("AR", n), df_ATD_temp)
        elseif num in enumerate(chem_cycle_list)
            for cyc in enumerate(chem_cycle_list)
                if num in cyc

                end                    
            end
        end
        
        df_ATD_temp = CES_checker(num, mol, mol_graph, df_ATD_temp, ATD_array)
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_array[num] = df_ATD_temp.type_name[1]
    end
    return ATD_array
end


function CES_checker(num::Integer, mol::Molecule, graph::SimpleGraph, df::DataFrame, ATD_array::AbstractArray)
    
    ### To Do: check for neighbors with Graphs. neighbors(graph, num) in 
    return df
end


function cycle_checker(graph::SimpleGraph)
    if !(ne(graph) >= nv(graph))
        println("no cycle detected")
        return false
    else
        # Graphs.jl cycle_basis function for all cycles
        all_cycles_list = cycle_basis(graph)
        ###

        return all_cycles_list
    end
end


function build_graph(mol::Molecule)
    mol_graph = SimpleGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i])
    end
    return mol_graph
end

