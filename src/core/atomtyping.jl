
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs

export select_atomtyping, get_atomtype, build_graph, cycle_checker, build_weighted_graph, toString


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


function get_atomtype(mol::AbstractMolecule, df_ATD::DataFrame)
    # Graph representation of Molecule
    mol_graph = build_graph(mol)
    adj_matrix = adjacency_matrix(mol_graph)
    mol_wgraph = build_weighted_graph(mol)
    wgraph_adj_matrix = adjacency_matrix(mol_wgraph)
       
    # Cycle detection and list
    cycle_bool = (ne(mol_graph) >= nv(mol_graph))
    chem_cycle_list = cycle_checker(mol_graph)
    ring_intersections_matrix = ring_Intersections(chem_cycle_list)
    ring_class_list = hueckel_test(chem_cycle_list)

    # Dataframe of Elements from Molecule that is later filled 
    # with specific atomtype definitions and returned
    # here: filling of Element_wNeighborCount and BondTypes
    ATD_df = DataFrame([Array{String, 1}(undef, nrow(mol.atoms)), Vector{Vector{String}}(undef, nrow(mol.atoms)), Array{String, 1}(undef, nrow(mol.atoms))], 
    ["Element_wNeighborCount", "BondTypes", "Possible_Atomtypes"])
    for i = (1:nrow(ATD_df)) 
        ATD_df.Element_wNeighborCount[i] = string(toString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        for (j, neigh) in enumerate(neighbors(mol_graph, i))
            bond_type_num = Int(wgraph_adj_matrix[i,neigh])
            if j == 1
                ATD_df.BondTypes[i] = toString(BondDef(bond_type_num))
            else
                ATD_df.BondTypes[i] = string(ATD_df.BondTypes[i], toString(BondDef(bond_type_num)))
                #push!(ATD_df.BondTypes[i], uppercase(toString(BondDef(bond_type_num))))
            end
        end
        for cycnum = (1:lastindex(chem_cycle_list)) 
            if cycle_bool && i in chem_cycle_list[cycnum]
                ATD_df.BondTypes[i] = lowercase(ATD_df.BondTypes[i])
            end 
        end
        ATD_df.BondTypes[i] = format_BondTypes(ATD_df.BondTypes[i])
    end

    for num = (1:nrow(mol.atoms))
        # filter for element in df
        element_num = element_number(num, mol)
        
        df_ATD_temp = filter(:atomic_number => x -> x == element_num, df_ATD)
        df_neighbors = neighbors(mol_graph, num)

        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
            continue
        end
        
        # Filtering if current element is H atom
        if element_num == 1
            neighbor_element = toString(mol.atoms.element[df_neighbors[1]])
            df_ATD_temp_save = copy(df_ATD_temp)
            df_ATD_temp = filter(:CES => n -> occursin(ATD_df.Element_wNeighborCount[neighbors(mol_graph, num)[1]], n), df_ATD_temp)
            if nrow(df_ATD_temp) == 0
                df_ATD_temp = copy(df_ATD_temp_save)
            end
            df_ATD_temp = filter(:CES => n -> occursin(neighbor_element, n), df_ATD_temp)
            elec_wgroups = count_withdrawal_groups(neighbors(mol_graph, num)[1], mol, mol_graph)
            df_ATD_temp = filter(:electron_withdrawal_groups => n -> n == elec_wgroups, df_ATD_temp)
            if neighbor_element == "O"
                ### check weither H20 classify as hw or not then classify as ho  
            elseif nrow(df_ATD_temp) > 1 
                
            end
        end
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
            continue
        end

        # filter for connected_atoms and connected_H_atoms
        connected_atoms = 0
        connected_H_atoms = 0
        for i = (1:nrow(mol.atoms))      
            if adj_matrix[num,i] >= 1
                connected_atoms += 1
                if element_number(i,mol) == 1
                    connected_H_atoms += 1
                end
            end
        end
        
        df_ATD_temp = filter(:num_neighbors => n -> n == connected_atoms, df_ATD_temp)
        if connected_H_atoms >= 1 && (0 in df_ATD_temp.num_H_bonds || 1 in df_ATD_temp.num_H_bonds || 2 in df_ATD_temp.num_H_bonds || 3 in df_ATD_temp.num_H_bonds)
            df_ATD_temp = filter(:num_H_bonds => n -> n == connected_H_atoms, df_ATD_temp)
        end

        # filter out obvious loop properties if no cycle detected
        # or add cycle length property through RG3,...,9 
        part_of_num_rings = 0
        if !occursin("AR", ATD_df.BondTypes[num]) && !occursin("RG", ATD_df.BondTypes[num])
            df_ATD_temp = filter(:atomic_property => n -> !occursin("RG", n), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> !occursin("AR", n), df_ATD_temp)
        end
        
        
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
        end
    end
    return ATD_df
end


function ring_Intersections(LList::Vector{Vector{Int64}})
    inters_matrix = Matrix{Vector{Int64}}(undef, lastindex(LList), lastindex(LList))
    for i = (1:lastindex(LList))
        curr1_list = LList[i]
        for j = (1:lastindex(LList))
            curr2_list = LList[j]
            inters_matrix[i,j] = inters_matrix[j,i] = intersect!(curr1_list, curr2_list)
        end
    end

function hueckel_test(LList::Vector{Vector{Int64}}, wgraph_adj::adjacency_matrix)
    ### n = (num_pi_elec - 2) / 4 , if n % 1 == 0 then Hueckel
    ring_list = Vector{Vector{String}}()
    for vlist in LList
        push!(ring_list, string("RG", String(lastindex(vlist))))
        pi_elec = 0
        for bond = (1:lastindex(vlist)-1)
            if bond == 1 && Int(wgraph_adj[vlist[bond], vlist[lastindex(vlist)]]) == 2
                pi_elec += 2
            elseif Int(wgraph_adj[vlist[bond], vlist[bond+1]]) == 2
                pi_elec += 2
            end
            bondsum += wgraph_adj_matrix[ringlist[bond], ringlist[lastindex(ringlist)]]
        end 
            for (i,vert) in enumerate(vlist)
                if num in ringlist
                    part_of_num_rings += 1
                    bondsum = 0
                    ring_size = lastindex(ringlist)
                    
                        
            end
        end
    end
    return
end


function count_withdrawal_groups(num::Int, mol::AbstractMolecule, mol_graph::SimpleGraph)
    elec_pullers = ["COOH", "OH", "SOO", "Cl", "F", "Br", "I", "NOO"] ### ToDo: extend for all EWD Groups
    elec_pullers_num = 0
    for neigh1 in neighbors(mol_graph, num)
        neighbor_str = toString(mol.atoms.element[neigh1])
        for neigh2 in neighbors(mol_graph, neigh1)
            neighbor_str = string(neighbor_str, toString(mol.atoms.element[neigh2]))
        end
        if neighbor_str in enumerate(elec_pullers)
            elec_pullers_num += 1
        end
    end
    if elec_pullers_num != 0
        return elec_pullers_num
    end
    return Int(-1)
end


function CES_checker(num::Integer, mol::AbstractMolecule, graph::SimpleGraph, df::DataFrame, ATD_array::AbstractArray)
    
    ### To Do: check for neighbors with Graphs. neighbors(graph, num) in 
    return df
end


function cycle_checker(graph::SimpleGraph)
    if !(ne(graph) >= nv(graph))
        println("no cycle detected")
        return [[]]
    else
        # Graphs.jl cycle_basis function for all cycles
        all_cycles_list = cycle_basis(graph)
        ###

        return all_cycles_list
    end
end


function build_graph(mol::AbstractMolecule)
    mol_graph = SimpleGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i])
    end
    return mol_graph
end


function build_weighted_graph(mol::AbstractMolecule)
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


function element_number(num::Int64, mol::AbstractMolecule)
    ele_num = Int8(mol.atoms.element[num])
    return ele_num
end


function toString(AnyEnum::Enum)
    return String(Symbol(AnyEnum))
end