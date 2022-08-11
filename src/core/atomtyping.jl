
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

    # Dataframe of Elements from Molecule that is later filled 
    # with specific atomtype definitions and returned
    # here: filling of Element_wNeighborCount and BondTypes
    ATD_df = DataFrame([Array{String, 1}(undef, nrow(mol.atoms)) for _ = (1:3)], ["Element_wNeighborCount", "BondTypes", "Possible_Atomtypes"])
    for i = (1:nrow(ATD_df)) 
        ATD_df.Element_wNeighborCount[i] = string(toString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        for (j, neigh) in enumerate(neighbors(mol_graph, i))
            bond_type_num = Int(wgraph_adj_matrix[i,neigh])
            if j == 1
                ATD_df.BondTypes[i] = toString(BondDef(bond_type_num))
            else
                ATD_df.BondTypes[i] = string(ATD_df.BondTypes[i], toString(BondDef(bond_type_num)))
            end    
        end
        if !(i in enumerate(chem_cycle_list))
            uppercase(ATD_df.BondTypes[i])
        end
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
        if !(num in enumerate(chem_cycle_list)) 
            df_ATD_temp = filter(:atomic_property => n -> !occursin("RG", n), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> !occursin("AR", n), df_ATD_temp)
        elseif num in enumerate(chem_cycle_list)
            for ringlist in enumerate(chem_cycle_list)
                if num in ringlist
                    part_of_num_rings += 1
                    bondsum = 0
                    ring_size = lastindex(ringlist)
                    for bond = (1:lastindex(ringlist)-1)
                        if bond == 1
                            bondsum += wgraph_adj_matrix[ringlist[bond], ringlist[lastindex(ringlist)]]
                        end
                        bondsum += wgraph_adj_matrix[ringlist[bond], ringlist[lastindex(ringlist)]]
                    end
                end                    
            end
        end
        
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
        end
    end
    return ATD_df
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