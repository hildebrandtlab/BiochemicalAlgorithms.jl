
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs

export select_atomtyping, get_atomtype, build_graph, build_graph_unweighted

# export atomtype_gaff

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
    
    # Cycle detection and number
    cycle_bool = (ne(mol_graph) >= nv(mol_graph))
    chem_cycle_list = cycle_checker(mol_graph)
    # cycle_number = (ne(mol_graph) - nv(mol_graph))+1


    for num = (1:nrow(mol.atoms))
        # filter for element in df
        element_gaff = Int8((filter(:number => n -> n == num, mol.atoms).element[1]))
        df_ATD_temp = filter(:atomic_number => x -> x == element_gaff, df_ATD)
        
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
        if !cycle_bool # || (!(2 in adjacency_matrix(graph)) && !(4 in adjacency_matrix(graph)))
            df_ATD_temp = filter(:atomic_property => n -> !occursin("RG", n), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> !occursin("AR", n), df_ATD_temp)
        end
        
        df_ATD_temp = CES_checker(num, mol, mol_graph, df_ATD_temp, ATD_array)
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_array[num] = df_ATD_temp.type_name[1]
        end   
    end
    return ATD_array
end


function cycle_checker(graph::SimpleWeightedGraph)
    if !(ne(graph) >= nv(graph))
        println("no cycle detected")
        return false
    else
        df_list = Array{DataFrame}

        ### Petitjean et al method for number of cycles
        cycle_count = ne(graph)-nv(graph)
        vertex_adj_count = zeros(Int64, 6)
        for i = (1:nv(graph))
            vdegree = degree(graph, i)
            vertex_adj_count[vdegree] += 1
        end
        println(vertex_adj_count)
        for j = (3:lastindex(vertex_adj_count))
            cycle_count += (((j-2)/2) * vertex_adj_count[j])
        end
 
        print(cycle_count)
        testing = readline()
        ###

        ### vertex removal method
        cycle_vertices_list = Vector{Int64}()
        for i = (1:nv(graph))
            temp_graph = graph
            rem_vertex!(temp_graph,i)
            # results in broken graph that sadly doesn't help finding cycles
            if (ne(temp_graph)-nv(temp_graph)) <= ((ne(graph)-nv(graph))+1)
                append!(cycle_vertices_list, i)
            end
        end
        ###
        
        ### spanning tree method
        missing_edge_list = Vector{}()
        # Fileformat problem with SimpleWeightedEdge
        mst_graph = boruvka_mst(graph)
        for i in edges(graph)
            if i in mst_graph[1]
                append!(missing_edge_list, i) 
            end
        end
        println(missing_edge_list)
        testing = readline()
        ###

        x_temp = ne(graph) - nv(graph) + 1
        println("naive: graph has at least",x_temp, " cycle(s)")
        testing = readline()
        return true
        # To Do: return list of cycles
    end
end

function CES_checker(num::Integer, mol::Molecule, graph::SimpleWeightedGraph, df::DataFrame, ATD_array::AbstractArray)
    
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
        elseif Int8(mol.bonds.order[i]) == 5
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 5)
        elseif Int8(mol.bonds.order[i]) == 6
            add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i], 6)
        end
    end
    return mol_graph
end

function build_graph_unweighted(mol::Molecule)
    mol_graph = SimpleGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i])
    end
    return mol_graph
end

