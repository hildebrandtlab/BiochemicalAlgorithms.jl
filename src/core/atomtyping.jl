
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs
using StatsBase

export select_atomtyping, get_atomtype, build_graph, cycle_checker, build_weighted_graph, enumToString, cycle_intersections


function select_atomtyping()

    ### temporary for faster testing
    df_ATD = CSV.read("atomtype_gff.csv", DataFrame)
    return df_ATD
    ###

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
    ring_class_list = NG_RG_AR_DEFtype(chem_cycle_list, wgraph_adj_matrix, ring_intersections_matrix, mol)
    
    # Dataframe of Elements from Molecule that is later filled 
    # with specific atomtype definitions and returned
    # here: filling of Element_wNeighborCount and BondTypes
    ATD_df = DataFrame([Array{String, 1}(undef, nrow(mol.atoms)), 
                Vector{Vector{String}}(undef, nrow(mol.atoms)), Vector{Vector{Int64}}(undef, nrow(mol.atoms)), 
                Vector{Vector{Vector{Int64}}}(undef, nrow(mol.atoms)), Vector{Vector{String}}(undef, nrow(mol.atoms))], 
                ["Element_wNeighborCount", "BondTypes", "Neighbors", "Secondary_Neighbors" ,"Possible_Atomtypes"])
    for i = (1:nrow(ATD_df)) 
        str_for_BondTypes = ""
        ATD_df.Element_wNeighborCount[i] = string(enumToString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        for neigh in neighbors(mol_graph, i)
            bond_type_num = Int(wgraph_adj_matrix[i,neigh])
            str_for_BondTypes = string(str_for_BondTypes, enumToString(DefBond(bond_type_num)))
        end
        ATD_df.BondTypes[i] = format_BondTypes!(i, str_for_BondTypes, ring_class_list)
        # neighbors and neighbors of neighbors for each atom into ATD_df.Neighbors
        ATD_df.Neighbors[i] = neighbors(mol_graph, i)
        ATD_df.Secondary_Neighbors[i] = Vector{Vector{Int64}}()
        for (n, prim_neigh) in enumerate(ATD_df.Neighbors[i])
            push!(ATD_df.Secondary_Neighbors[i],[]) #Careful, neighbor of neighbor of H, Halogens, double Bond O and S are now empty lists
            for sec_neigh in neighbors(mol_graph, prim_neigh)
                if sec_neigh != i
                    push!(ATD_df.Secondary_Neighbors[i][n], sec_neigh)
                end
            end
        end
    end
    
    # Filter Dataframe process for each atom in molecule
    for num = (1:nrow(mol.atoms))
        # filter for element in df
        element_num = element_number(num, mol)
        df_ATD_temp = filter(:atomic_number => x -> x == element_num, df_ATD)

        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
            continue
        end
        
        # Filtering if current element is H atom
        if element_num == 1
            df_ATD_temp_save = copy(df_ATD_temp)
            # filter for electron withdrawing groups and Neighbor-element
            EWG_num = count_EWG(num, ATD_df)
            df_ATD_temp = filter(:electron_withdrawal_groups => n -> n == EWG_num, df_ATD_temp)
            neighbor_element = ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num][1]]
            df_ATD_temp = filter(:CES => n -> occursin(neighbor_element[1], n), df_ATD_temp)
            if nrow(df_ATD_temp) == 1
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            elseif nrow(df_ATD_temp) == 0
                df_ATD_temp = copy(df_ATD_temp_save)
            end
            
            # case: neighbor is Oxygen
            if neighbor_element == "O2" && in(ATD_df.Secondary_Neighbors[num][1]).("H1")
                df_ATD_temp = filter(:CES => n -> n == "(O(H1))", df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            elseif neighbor_element == "O2" && !in(ATD_df.Secondary_Neighbors[num][1]).("H1")
                df_ATD_temp = filter(:CES => n -> n == "(O)", df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            end
            
            # case: neighbor is Carbon
            if occursin("C", neighbor_element) && !in(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][1]]).("N4")
                df_ATD_temp = filter(:CES => n -> occursin(neighbor_element, n), df_ATD_temp)
            elseif occursin("C", neighbor_element) && in(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num]]).("N4")
                df_ATD_temp = filter(:CES => n -> n == "(C(N4))", df_ATD_temp)
            end
            if nrow(df_ATD_temp) == 1
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            elseif nrow(df_ATD_temp) == 0
                df_ATD_temp = copy(df_ATD_temp_save)
            end

            # other cases
            if occursin("N", neighbor_element)
                df_ATD_temp = filter(:CES => n -> n == "(N)" , df_ATD_temp)
            elseif occursin("S", neighbor_element)
                df_ATD_temp = filter(:CES => n -> n == "(S)" , df_ATD_temp)
            elseif occursin("P", neighbor_element)
                df_ATD_temp = filter(:CES => n -> n == "(P)" , df_ATD_temp)
            else
                df_ATD_temp = filter(:CES => n -> (occursin(neighbor_element, n) || occursin("X", n)) , df_ATD_temp)
            end
            if nrow(df_ATD_temp) == 1
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            elseif nrow(df_ATD_temp) == 0
                df_ATD_temp = filter(:type_name => n -> n == "ha" , df_ATD_temp_save)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            end
        end

        # filter for connected_atoms
        connected_atoms = lastindex(neighbors(mol_graph, num))
        df_ATD_temp = filter(:num_neighbors => n -> n == connected_atoms, df_ATD_temp)
        
        # filter for connected_H_atoms
        if (true in in(df_ATD_temp.num_H_bonds).([1,2,3]) && in(mol.atoms.element[ATD_df.Neighbors[num]]).(Elements.H)) ||
            (in(df_ATD_temp.num_H_bonds).(0) && !in(mol.atoms.element[ATD_df.Neighbors[num]]).(Elements.H))
            connected_H_atoms = countmap(in([Elements.H]).(mol.atoms.element[ATD_df.Neighbors[num]]))[true]
            df_ATD_temp = filter(:num_H_bonds => n -> n == connected_H_atoms, df_ATD_temp)
        end

        # filter out obvious loop properties if no cycle detected
        # To Do: filter for all loop properties in BondTypes, maybe finding and selecting through properties instead of process of elimination as seen below
        AR1_RG6 = ["RG6", "AR1"]
        AR2_RG6 = ["RG6", "AR2"]
        if "NG" in ATD_df.BondTypes[num]
            df_ATD_temp = filter(:atomic_property => n -> (!occursin("RG", n) && !occursin("AR", n)), df_ATD_temp)
            df_ATD_temp = filter(:atomic_property => n -> (!occursin("sb", n) && !occursin("db", n)), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR1") #&& in(ATD_df.BondTypes[neigh_list]).("NG")
            # nur NG rausfiltern funktioniert nicht für GAFF def. Und in BCC DEF nicht eindeutig (s. Zeile 18)
            df_ATD_temp = filter(:atomic_property => n -> (occursin("AR1", n) && !occursin("RG6", n)), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR2")
            df_ATD_temp = filter(:atomic_property => n -> occursin("AR2", n), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR3")
            df_ATD_temp = filter(:atomic_property => n -> occursin("AR3", n), df_ATD_temp)
        end

        # specific non-cycle C3 filter for XA1 (which is O1 or S1) or XB (which are the Halogens)
        XB_elements = [Elements.Cl, Elements.Br, Elements.F, Elements.I, Elements.At, Elements.Ts]
        XA1_elements = [Elements.O, Elements.S]
        XB_ATD_df = ["Cl1", "F1", "I1", "Br1", "At1", "Ts1"]
        XA1_ATD_df = ["O1", "S1"]
        df_ATD_temp_save = copy(df_ATD_temp)
        if ATD_df.Element_wNeighborCount[num] == "C3" && "NG" in ATD_df.BondTypes[num]
            if true in in(mol.atoms.element[ATD_df.Neighbors[num]]).(XB_elements)  
                df_ATD_temp = filter(:CES => n -> (occursin("XB", n) || occursin("*", n)), df_ATD_temp)
            elseif true in in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).(XA1_ATD_df)
                df_ATD_temp = filter(:CES => n -> occursin("XA1", n), df_ATD_temp)
            elseif countmap(in(["N3"]).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]))[1] == 3
                df_ATD_temp = filter(:CES => n -> n == "(N3,N3,N3)", df_ATD_temp)
            end
        end
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
            continue
        elseif nrow(df_ATD_temp) == 0
            df_ATD_temp = copy(df_ATD_temp_save)
        end

        # specific non-cycle N3 filtering
        if ATD_df.Element_wNeighborCount[num] == "N3" && "NG" in ATD_df.BondTypes[num]
            if countmap(in(["C3"]).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]))[1] == 1 

            elseif in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).("C4") && !in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(2)
                df_ATD_temp = filter(:CES => n -> n == "(C3(XA1))", df_ATD_temp)
            elseif true in in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).(XA1_ATD_df)
                df_ATD_temp = filter(:CES => n -> occursin("XA1", n), df_ATD_temp)
            elseif countmap(in(["N3"]).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]))[1] == 3
                df_ATD_temp = filter(:CES => n -> n == "(N3,N3,N3)", df_ATD_temp)
            end
        end
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
            continue
        elseif nrow(df_ATD_temp) == 0
            df_ATD_temp = copy(df_ATD_temp_save)
        end

        # add all left over atomtypes into a list
        if nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
            continue
        else
            count_dict = countmap(df_ATD_temp.type_name)
            ATD_df.Possible_Atomtypes[num] = []
            for (key, value) in count_dict
                push!(ATD_df.Possible_Atomtypes[num], key)
            end
        end
    end
    return ATD_df
end


function c_tag(num::Int64, ATD_df::DataFrame)
    # the "c" tag is added in GAFF.DEF to atomtypes in aromatic rings with one of the following properties
    bool_c_tag = false
    c_tag_list = [["(C3(C3))"], ["(C3(C2))"], ["(C3(XB3))"], ["(XB2(XB2))"], ["(XB2(C2))"], ["(XB2(C3))"],
                    ["(C3[sb'])"], ["(XB2[sb'])"], ["(XD3[sb',db])"], ["(XD4[sb',db])"]]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        for sec_neigh in ATD_df.Secondary_Neighbors[num][i]

            if in() ### build neighbor string and check if in c_tag_list

            end
        end
    end
    return bool_c_tag
end


function format_BondTypes!(num::Int64, str_Bonds::AbstractString, ring_class_list::Vector{Vector{String}})
    str_bondtypes_list = Vector{String}()
    non_ring_atom_bool = ("NG" in ring_class_list[num]) 
    for defi = (1:3)
        push!(str_bondtypes_list, enumToString(DefBond(defi)))
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


function NG_RG_AR_DEFtype(LList::Vector{Vector{Int64}}, wgraph_adj::Graphs.SparseMatrix, inters_matrix::Matrix{Vector{Int64}}, mol::AbstractMolecule)
    ### n = (num_db*2 - 2) / 4 , if n % 1 == 0 then Hueckel
    ring_class_list = Vector{Vector{String}}(undef, nrow(mol.atoms))
    for i = (1:nrow(mol.atoms))
        ring_class_list[i] = ["NG"]
    end
    for (numvlist, vlist) in enumerate(LList)
        for x in vlist
            if ring_class_list[x] == ["NG"]
                ring_class_list[x] = [string("RG", string(lastindex(vlist)))]
            end
        end

        # check if is O, N, or S present in Ring vertex list
        ONS_present = false
        if true in in(mol.atoms.element[vlist]).([Elements.O,Elements.N,Elements.S])
            ONS_present = true
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
        if (pi_elec / lastindex(vlist)) == 1.0 && !ONS_present
            for x in vlist
                push!(ring_class_list[x], "AR1")
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) <= 1 && ONS_present && lastindex(vlist) > 4
            for x in vlist
                push!(ring_class_list[x], "AR2")
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) < 1 && !ONS_present && lastindex(vlist) > 4
            # check if Ring vlist has intersections with other rings in molecule and if these are aromatic
            #intersection_db = countmap(in(inters_matrix).()
            has_aromatic_inters = false
            for i = (1:lastindex(LList))
                if !isempty(inters_matrix[numvlist,i]) && lastindex(inters_matrix[numvlist, i]) == 2
                    atom1_bonds = filter(:a1 => n -> n == inters_matrix[numvlist,i][1], mol.bonds)
                    atom2_bonds = filter(:a1 => n -> n == inters_matrix[numvlist,i][2], mol.bonds)
                    if countmap(in([BondOrder.Double]).(atom1_bonds.order))[1] == 1 &&
                        countmap(in([BondOrder.Double]).(atom2_bonds.order))[1] == 1
                        for x in vlist
                            if !in(ring_class_list[x]).("AR1")
                                push!(ring_class_list[x], "AR1")
                            end
                        end 
                    end           
                end
            end
        elseif (pi_elec / lastindex(vlist)) == 0
            for x in vlist
                if !in(ring_class_list[x]).("AR5")
                    push!(ring_class_list[x], "AR5")    
                end
            end
        end
    end
    return ring_class_list
end


function count_EWG(num::Int64, ATD_df::DataFrame)
    strong_pullers = ["Cl1", "F1", "Br1", "I1", "O1", "S1"]
    possible_pullers = ["C2", "C3", "C4", "S3", "N3", "P3", "P4", "O2", "S2"] ### To Do: extend for all EWD Groups
    elec_pullers_num = 0
    if true in in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).(possible_pullers)
        for i = (1:lastindex(ATD_df.Neighbors[num]))
            if true in in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]])
                # all neighbors of neighbor that are in strong_pullers are an EWG
                elec_pullers_num += countmap(in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]]))[true]
            end
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if in(possible_pullers).(ATD_df.Element_wNeighborCount[sec_neigh]) && 
                    true in in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[sec_neigh]])
                    # For example for NO2 group. Only counts as one EWG
                    elec_pullers_num += 1
                end  
            end  
        end
    end
    if elec_pullers_num == 0
        elec_pullers_num = -1
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


function enumToString(AnyEnum::Enum)
    return String(Symbol(AnyEnum))
end