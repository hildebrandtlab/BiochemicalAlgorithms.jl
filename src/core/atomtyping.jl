
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs

export select_atomtyping, get_atomtype, build_graph, cycle_checker, build_weighted_graph, toString, cycle_intersections


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
    ring_intersections_matrix = Matrix{Vector{Int64}}(undef, lastindex(chem_cycle_list), lastindex(chem_cycle_list))
    ring_class_list = Vector{Vector{String}}()
    ring_intersections_matrix = cycle_intersections(chem_cycle_list)
    ring_class_list = NG_RG_AR_DEFtype(chem_cycle_list, mol_graph, wgraph_adj_matrix, ring_intersections_matrix, mol)
    
    
    # Dataframe of Elements from Molecule that is later filled 
    # with specific atomtype definitions and returned
    # here: filling of Element_wNeighborCount and BondTypes
    ATD_df = DataFrame([Array{String, 1}(undef, nrow(mol.atoms)), Vector{Vector{String}}(undef, nrow(mol.atoms)), Array{String, 1}(undef, nrow(mol.atoms))], 
    ["Element_wNeighborCount", "BondTypes", "Possible_Atomtypes"])
    for i = (1:nrow(ATD_df)) 
        str_for_BondTypes = ""
        ATD_df.Element_wNeighborCount[i] = string(toString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        for neigh in neighbors(mol_graph, i)
            bond_type_num = Int(wgraph_adj_matrix[i,neigh])
            str_for_BondTypes = string(str_for_BondTypes, toString(BondDef(bond_type_num)))
        end
        ATD_df.BondTypes[i] = format_BondTypes!(i, str_for_BondTypes, ring_class_list)
    end
    
    # Filter Dataframe process for each atom in molecule
    for num = (1:nrow(mol.atoms))
        # filter for element in df
        element_num = element_number(num, mol)
        df_ATD_temp = filter(:atomic_number => x -> x == element_num, df_ATD)
        neigh_list = neighbors(mol_graph, num)

        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
            continue
        end
        
        # Filtering if current element is H atom
        if element_num == 1
            neighbor_element = toString(mol.atoms.element[neigh_list[1]])
            elec_wgroups = count_withdrawal_groups(neighbors(mol_graph, num)[1], mol, mol_graph)
            df_ATD_temp_save = copy(df_ATD_temp)
            # case: neighbor is Oxygen
            if neighbor_element == "O" && all(in(["H1"]).(ATD_df.Element_wNeighborCount[neighbors(mol_graph,neigh_list)]))
                df_ATD_temp = filter(:CES => n -> n == "(O(H1))", df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            elseif neighbor_element == "O" && !all(in(["H1"]).(ATD_df.Element_wNeighborCount[neighbors(mol_graph,neigh_list)]))
                df_ATD_temp = filter(:CES => n -> n == "(O)", df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            end
            # case: C3 with no EWD
            if ATD_df.Element_wNeighborCount[neigh_list[1]] == "C3" && elec_wgroups == -1
                df_ATD_temp = filter(:CES => n -> n == "*", df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            elseif ATD_df.Element_wNeighborCount[neigh_list[1]] == "C3" && elec_wgroups >= 1
                df_ATD_temp = filter(:CES => n -> (occursin(ATD_df.Element_wNeighborCount[neigh_list[1]], n) || occursin("XX", n)), df_ATD_temp)
                df_ATD_temp = filter(:electron_withdrawal_groups => n -> n == elec_wgroups, df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            end

            # other cases
            df_ATD_temp = filter(:CES => n -> (occursin(neighbor_element, n) || occursin("X", n)) , df_ATD_temp)
            if nrow(df_ATD_temp) == 1
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            elseif nrow(df_ATD_temp) == 0
                df_ATD_temp = copy(df_ATD_temp_save)
            end

            df_ATD_temp = filter(:CES => n -> (occursin(ATD_df.Element_wNeighborCount[neigh_list[1]], n) || occursin("XX", n)), df_ATD_temp)
            df_ATD_temp = filter(:electron_withdrawal_groups => n -> n == elec_wgroups, df_ATD_temp)
            if nrow(df_ATD_temp) == 1
                ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
                continue
            elseif nrow(df_ATD_temp) == 0
                df_ATD_temp = copy(df_ATD_temp_save)
            end
        end

        # filter for connected_atoms
        connected_atoms = lastindex(neighbors(mol_graph, num))
        df_ATD_temp = filter(:num_neighbors => n -> n == connected_atoms, df_ATD_temp)
        
        # filter for connected_H_atoms
        if true in in(df_ATD_temp.num_H_bonds).([1,2,3]) && in(mol.atoms.element[neigh_list]).(Elements.H)
            connected_H_atoms = countmap(mol.atoms.element[neigh_list])[Elements.H]
            df_ATD_temp = filter(:num_H_bonds => n -> n == connected_H_atoms, df_ATD_temp)
        end

        # filter out obvious loop properties if no cycle detected
        # To Do: filter for all loop properties in BondTypes, maybe finding and selecting through properties instead of process of elimination as seen below
        if "NG" in ATD_df.BondTypes[num]
            df_ATD_temp = filter(:atomic_property => n -> (!occursin("RG", n) && !occursin("AR", n)), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> (!occursin("sb", n) && !occursin("db", n)), df_ATD_temp)
        end
        
        println(df_ATD_temp)
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = df_ATD_temp.type_name[1]
        end
    end
    return ATD_df
end


function format_BondTypes!(num::Int64, str_Bonds::AbstractString, ring_class_list::Vector{Vector{String}})
    str_bondtypes_list = Vector{String}()
    non_ring_atom_bool = ("NG" in ring_class_list[num]) 
    for defi = (1:3)
        push!(str_bondtypes_list, toString(BondDef(defi)))
    end
    ret_list = Vector{String}()
    for (j, type) in enumerate(str_bondtypes_list)
        count_type = count(==(type[1]), str_Bonds)
        if count_type > 0 && !non_ring_atom_bool
            push!(ret_list, string(count_type, type))  
        elseif count_type > 0 && non_ring_atom_bool
            push!(ret_list, string(count_type, uppercase(type)))
        end
    end
    append!(ret_list, ring_class_list[num])
    return ret_list
end


function NG_RG_AR_DEFtype(LList::Vector{Vector{Int64}}, mol_graph::SimpleGraph, wgraph_adj::Graphs.SparseMatrix, inters_matrix::Matrix{Vector{Int64}}, mol::AbstractMolecule)
    ### n = (num_db*2 - 2) / 4 , if n % 1 == 0 then Hueckel
    ring_class_list = Vector{Vector{String}}(undef, nrow(mol.atoms))
    for i = (1:nrow(mol.atoms))
        ring_class_list[i] = ["NG"]
    end
    for (numvlist, vlist) in enumerate(LList)
        for x in vlist
            if ring_class_list[x] == ["NG"]
                ring_class_list[x] = ["RG"]
                push!(ring_class_list[x], "AR")
                push!(ring_class_list[x], string("RG", string(lastindex(vlist))))
            end
        end

        # check if is O, N, or S present in Ring vlist
        ONS_present = false
        for i = (1:lastindex(vlist))
            if Int(mol.atoms.element[vlist[i]]) == 8 || Int(mol.atoms.element[vlist[i]]) == 7 || Int(mol.atoms.element[vlist[i]]) == 16
                ONS_present = true
            end
        end

        # check number of pi electrons
        pi_elec = 0
        for bond = (1:lastindex(vlist)-1)
            if bond == 1 && Int(wgraph_adj[vlist[bond], vlist[lastindex(vlist)]]) == 2
                pi_elec += 2
            elseif Int(wgraph_adj[vlist[bond], vlist[bond+1]]) == 2
                pi_elec += 2
            end
        end
        if (pi_elec / lastindex(vlist)) == 1.0
            for x in vlist
                push!(ring_class_list[x], "AR1")
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) < 1 && ONS_present
            for x in vlist
                push!(ring_class_list[x], "AR2")
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) < 1 && !ONS_present
            # check if Ring vlist has intersections with other rings in molecule and if these are aromatic
            has_aromatic_inters = false
            for i = (1:lastindex(LList))
                if !isempty(inters_matrix[numvlist,i]) && lastindex(inters_matrix[numvlist, i]) == 2
                    for j in inters_matrix[numvlist,i]
                        db_inters_atom = 0
                        for k in neighbors(mol_graph, j)
                            if wgraph_adj[j,k] == 2
                                db_inters_atom += 1
                            end
                        end
                        if db_inters_atom == 1
                            has_aromatic_inters = true
                            push!(ring_class_list[j], "AR3")
                        end                    
                    end                 
                end
            end
        elseif (pi_elec / lastindex(vlist)) < 1/2
            for x in vlist
                push!(ring_class_list[x], "AR5")
            end
        end
    end
    return ring_class_list
end


function count_withdrawal_groups(num::Int, mol::AbstractMolecule, mol_graph::SimpleGraph)
    # rework to create suitable neighborstring or if cascade for certain elements
    elec_pullers = ["COO", "CHO", "COH", "SOO", "Cl", "F", "Br", "I", "NOO"] ### To Do: extend for all EWD Groups
    elec_pullers_num = 0
    for neigh1 in neighbors(mol_graph, num)
        neighbor_str = toString(mol.atoms.element[neigh1])
        for neigh2 in neighbors(mol_graph, neigh1)
            neighbor_str = string(neighbor_str, toString(mol.atoms.element[neigh2]))
        end
        if in(elec_pullers).(neighbor_str)
            println("we got an elec_puller")
            elec_pullers_num += 1
        end
    end
    if elec_pullers_num == 0
        return Int(-1)
    end
    return elec_pullers_num
end


function cycle_intersections(LList::Vector{Vector{Int64}})
    inters_matrix = Matrix{Vector{Int64}}(undef, lastindex(LList), lastindex(LList))
    for i = (1:lastindex(LList))
        curr1_list = copy(LList[i])
        for j = (1:lastindex(LList))
            curr2_list = copy(LList[j])
            if curr1_list != curr2_list
                inters_matrix[i,j] = inters_matrix[j,i] = intersect(curr1_list, curr2_list)
            else
                inters_matrix[i,j] = inters_matrix[j,i] = []
            end
        end
    end
    return inters_matrix
end


function cycle_checker(graph::SimpleGraph)
    if !(ne(graph) >= nv(graph))
        println("no cycle detected")
        empty_vec = Vector{Vector{Int64}}()
        return empty_vec
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