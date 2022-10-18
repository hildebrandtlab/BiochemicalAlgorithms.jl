
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs
using StatsBase

export select_atomtyping, get_atomtype, build_graph, cycle_checker, build_weighted_graph, enumToString, cycle_intersections


function select_atomtyping()

    ### temporary for faster testing
    df_ATD = CSV.read("atomtype_gff.csv", DataFrame)
    return df_ATD
    ###

    # println("Which atomtyping do you want? \n Press: 1 for GAFF,  2 for Amber")
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
    # here: filling of Element_wNeighborCount, BondTypes, Neighbors and Secondary_Neighbors
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

        if isassigned(ATD_df.Possible_Atomtypes, num)
            if lastindex(ATD_df.Possible_Atomtypes[num]) == 1
                continue
            end
        elseif nrow(df_ATD_temp) == 1
            ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
            continue
        end

        neighbor_element_list = Vector{String}()
        for neigh in neighbors(mol_graph, num)
            push!(neighbor_element_list, enumToString(mol.atoms.element[neigh]))
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
            elseif (nrow(df_ATD_temp) == 0 || (nrow(df_ATD_temp) == 2 && EWG_num == -1)) && ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num][1]] == "C3"
                df_ATD_temp = filter(:type_name => n -> n == "ha", df_ATD_temp_save)
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
            elseif occursin("C", neighbor_element) && in(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][1]]).("N4")
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
        connected_H_atoms = -1
        if (true in in(df_ATD_temp.num_H_bonds).([1,2,3]) && in(mol.atoms.element[ATD_df.Neighbors[num]]).(Elements.H)) ||
            (in(df_ATD_temp.num_H_bonds).(0) && !in(mol.atoms.element[ATD_df.Neighbors[num]]).(Elements.H))
            connected_H_atoms = countmap(in([Elements.H]).(mol.atoms.element[ATD_df.Neighbors[num]]))[true]
        end
        df_ATD_temp = filter(:num_H_bonds => n -> n == connected_H_atoms, df_ATD_temp)

        # specific C4 filtering
        if ATD_df.Element_wNeighborCount[num] == "C4"
            if in(ATD_df.BondTypes[num]).("RG3")
                df_ATD_temp = filter(:type_name => n -> n == "cx" , df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            elseif in(ATD_df.BondTypes[num]).("RG4")
                df_ATD_temp = filter(:type_name => n -> n == "cy" , df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            else
                df_ATD_temp = filter(:type_name => n -> n == "c3" , df_ATD_temp)
                ATD_df.Possible_Atomtypes[num] = [df_ATD_temp.type_name[1]]
                continue
            end
        end

        # tag functions for grouping of properties according to GAFF.DEF
        # wildcat Elements: X Types Dictionary, as seen in from antechamber documentation
        X_dict = Dict{String, Vector{String}}("XX"=>["C","N","O","S","P"], "XA"=>["O","S"], "XB"=>["N","P"], "XD"=>["S","P"])
        is_c_tagged, is_e_tagged, is_g_tagged, is_x_tagged, is_y_tagged, is_h_tagged = repeat([false],6)
        if ((in([Elements.N, Elements.P]).(mol.atoms.element[num]) && lastindex(ATD_df.Neighbors[num]) == 2) || 
                (in([Elements.C]).(mol.atoms.element[num]) && lastindex(ATD_df.Neighbors[num]) == 3))  && 
                (in(ATD_df.BondTypes[num]).("AR2") || in(ATD_df.BondTypes[num]).("AR3"))
            is_c_tagged = DEF_c_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", c_tag: ", is_c_tagged)
            if is_c_tagged
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "c")), df_ATD_temp)
            end
        elseif ((in([Elements.N, Elements.P]).(mol.atoms.element[num]) && lastindex(ATD_df.Neighbors[num]) == 2) || 
                (in([Elements.C]).(mol.atoms.element[num]) && lastindex(ATD_df.Neighbors[num]) == 3)) && 
                in(ATD_df.BondTypes[num]).("NG") && all(in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).([1,2]))
            is_e_tagged = DEF_e_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", e_tag: ", is_e_tagged)
            if is_e_tagged        
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "e")), df_ATD_temp)
            end
        elseif ATD_df.Element_wNeighborCount[num] == "C2" && all(in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).([1,3]))
            is_g_tagged = DEF_g_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", g_tag: ", is_g_tagged)
            if is_g_tagged
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "g")), df_ATD_temp)
            end
        elseif (ATD_df.Element_wNeighborCount[num] == "P3" || ATD_df.Element_wNeighborCount[num] == "S3") && in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(2)
            is_x_tagged = DEF_x_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", x_tag: ", is_x_tagged)
            if is_x_tagged
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "x")), df_ATD_temp)
            end
        elseif (ATD_df.Element_wNeighborCount[num] == "P4" || ATD_df.Element_wNeighborCount[num] == "S4") && in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(2)
            is_y_tagged = DEF_y_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", y_tag: ", is_y_tagged)
            if is_y_tagged
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "y")), df_ATD_temp)
            end
        elseif ATD_df.Element_wNeighborCount[num] == "N3" && true in in(X_dict["XX"]).(neighbor_element_list) && in(ATD_df.BondTypes[num]).("NG")
            is_h_tagged = DEF_h_tag(num, ATD_df, wgraph_adj_matrix, X_dict, mol)
            # println("num = ", num, ", h_tag: ", is_h_tagged)
            if is_h_tagged
                df_ATD_temp = filter(:type_name => n -> (lastindex(n) == 2 && n == string(lowercase(ATD_df.Element_wNeighborCount[num][1]), "h")), df_ATD_temp)
            end
        end
        if length(countmap(df_ATD_temp.type_name)) == 1
            count_dict = countmap(df_ATD_temp.type_name)
            ATD_df.Possible_Atomtypes[num] = []
            for (key, value) in count_dict
                push!(ATD_df.Possible_Atomtypes[num], key)
            end
            continue
        end

        # specific non-cycle C3 filter for XA1 (which is O1 or S1) or XB1 (which here are the Halogens)
        XB1_elements = [Elements.Cl, Elements.Br, Elements.F, Elements.I, Elements.At, Elements.Ts]
        XA1_elements = [Elements.O, Elements.S]
        XB_ATD_df = ["Cl1", "F1", "I1", "Br1", "At1", "Ts1"]
        XA1_ATD_df = ["O1", "S1"]
        df_ATD_temp_save = copy(df_ATD_temp)
        if ATD_df.Element_wNeighborCount[num] == "C3" && true in in(ATD_df.BondTypes[num]).(["NG", "AR5"]) && !is_e_tagged
            if true in in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).(XA1_ATD_df)
                df_ATD_temp = filter(:CES => n -> occursin("(XA1)", n), df_ATD_temp)
                # to get type_name "c" from GAFF.DEF
            elseif countmap(in(["N3"]).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]))[1] == 3
                df_ATD_temp = filter(:CES => n -> n == "(N3,N3,N3)", df_ATD_temp)
                # to get type_name "cz" from GAFF.DEF
            else
                df_ATD_temp = filter(:type_name => n -> n == "c2", df_ATD)
            end
        end

        # filter out obvious loop properties if no cycle detected
        # To Do: filter for all loop properties in BondTypes, maybe finding and selecting through properties instead of process of elimination as seen below
        if "NG" in ATD_df.BondTypes[num] && !is_e_tagged
            df_ATD_temp = filter(:atomic_property => n -> (!occursin("RG", n) && !occursin("AR", n) && !occursin("sb", n) && !occursin("db", n)), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR1") 
            df_ATD_temp = filter(:atomic_property => n -> occursin("AR1", n), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR2") && is_c_tagged
            df_ATD_temp = filter(:atomic_property => n -> occursin("AR2", n), df_ATD_temp)
        elseif in(ATD_df.BondTypes[num]).("AR3") && is_c_tagged
            df_ATD_temp = filter(:atomic_property => n -> occursin("AR3", n), df_ATD_temp)
        end

        # AR1 typed C3, cases: cp or ca
        if ATD_df.Element_wNeighborCount[num] == "C3" && in(ATD_df.BondTypes[num]).("AR1") &&
                all(in(["AR1"]).(ATD_df.BondTypes[ATD_df.Neighbors[num]]))
            df_ATD_temp = filter(:type_name => n -> n == "cp", df_ATD) 
        elseif ATD_df.Element_wNeighborCount[num] == "C3" && in(ATD_df.BondTypes[num]).("AR1") &&
            !all(in(["AR1"]).(ATD_df.BondTypes[ATD_df.Neighbors[num]]))
            df_ATD_temp = filter(:type_name => n -> n == "ca", df_ATD_temp)
        end

        # specific N3 and N2 filtering
        if ATD_df.Element_wNeighborCount[num] == "N3" && true in in(ATD_df.BondTypes[num]).(["NG", "AR5"]) 
            if in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).("C3") && countmap(in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).([2,3]))[1] == 0
                for (i,prim_neigh) in enumerate(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]])
                    if true in in(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]]).(XA1_ATD_df)
                        if in(ATD_df.BondTypes[num]).("RG3")
                            df_ATD_temp = filter(:atomic_property  => m -> m == "[RG3]", df_ATD_temp)
                            df_ATD_temp = filter(:CES => n -> n == "(C3(XA1))", df_ATD_temp)
                            # to get type_name "ni" from GAFF.DEF
                        elseif in(ATD_df.BondTypes[num]).("RG4")
                            df_ATD_temp = filter(:atomic_property  => m -> m == "[RG4]", df_ATD_temp)
                            df_ATD_temp = filter(:CES => n -> n == "(C3(XA1))", df_ATD_temp)
                            # to get type_name "nj" from GAFF.DEF
                        else
                            df_ATD_temp = filter(:type_name => n -> n == "n", df_ATD_temp)
                        end
                    end
                    # for sec_neigh in ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]]
                    #     println(sec_neigh)
                    #     if prim_neigh == "C3" && sec_neigh[lastindex(sec_neigh)] == "1"
                    #         df_ATD_temp = filter(:CES => n -> n == "(C3(XA1))", df_ATD_temp)
                    #         # to get type_name "n" from GAFF.DEF
                    #     end
                    # end
                end
            elseif countmap(in(["N3"]).(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]))[1] == 2
                df_ATD_temp = filter(:type_name => n -> n == "no", df_ATD)
            else
                df_ATD_temp = filter(:type_name => n -> n == "n3", df_ATD)
            end
        elseif ATD_df.Element_wNeighborCount[num] == "N2" && "NG" in ATD_df.BondTypes[num] && !(is_c_tagged || is_e_tagged)
            if all(in([2]).(wgraph_adj_matrix[num, ATD_df.Neighbors[num]])) || all(in([1,3]).(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]))
                df_ATD_temp = filter(:type_name => n -> n == "n1", df_ATD)
            elseif all(in([1,2]).(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]))
                df_ATD_temp = filter(:type_name => n -> n == "n2", df_ATD)
            end
        end

        # specific non-cycle S2, S3, and S4 filtering
        if ATD_df.Element_wNeighborCount[num] == "S2" 
            if (in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(2) || in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(3))
                df_ATD_temp = filter(:type_name => n -> n == "s2", df_ATD)
            elseif !in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(2) && !in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).(3)
                df_ATD_temp = filter(:type_name => n -> n == "ss", df_ATD)
            end
        elseif ATD_df.Element_wNeighborCount[num] == "S3" && !is_x_tagged
            df_ATD_temp = filter(:type_name => n -> n == "s4", df_ATD)
        elseif ATD_df.Element_wNeighborCount[num] == "S4" && !is_y_tagged
            df_ATD_temp = filter(:type_name => n -> n == "s6", df_ATD)
        end

        # ring based N3 and N2 filtering
        if ATD_df.Element_wNeighborCount[num] == "N3" && true in in(ATD_df.BondTypes[num]).(["AR1", "AR2", "AR3"]) && 
                countmap(in(wgraph_adj_matrix[num, ATD_df.Neighbors[num]]).([2,3]))[1] == 0 && !is_c_tagged
            df_ATD_temp = filter(:type_name => n -> n == "na", df_ATD_temp) 
        elseif ATD_df.Element_wNeighborCount[num] == "N2" && in(ATD_df.BondTypes[num]).("AR1")
            df_ATD_temp = filter(:type_name => n -> n == "nb", df_ATD_temp)
        end

        # Oxygen filtering, cases: op, oq, os
        if ATD_df.Element_wNeighborCount[num] == "O2" 
            if in(ATD_df.BondTypes[num]).("RG3")
                df_ATD_temp = filter(:type_name => n -> n == "op", df_ATD) 
            elseif in(ATD_df.BondTypes[num]).("RG4")
                df_ATD_temp = filter(:type_name => n -> n == "oq", df_ATD_temp)
            else
                df_ATD_temp = filter(:type_name => n -> n == "os", df_ATD_temp)
            end
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
    # println(ATD_df)
    return ATD_df
end


function DEF_c_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "c" tag is added in GAFF.DEF to atomtypes in aromatic rings with one or more of the c_tag_list properties
    c_tag_list = ["(C3(C3))", "(C3(C2))", "(C3(XB3))", "(XB2(XB2))", "(XB2(C2))", "(XB2(C3))",
                    "(C3[sb'])", "(XB2[sb'])", "(XD3[sb',db])", "(XD4[sb',db])"]
    c_tag_counter = 0 
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        if in(ATD_df.BondTypes[num]).("AR2") || in(ATD_df.BondTypes[num]).("AR3")
            # since newer version of GAFF.DEF c tag definition is deducible from preclassification of aromatic Ring type AR2 and AR3
            c_tag_counter += 2
            break
        elseif !in(ATD_df.BondTypes[prim_neigh]).("NG") && c_tag_counter < 2
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if !in(ATD_df.BondTypes[sec_neigh]).("NG")
                    # build neighbor-"neighbors of neighbor" or bond tuple and check if in c_tag_list
                    # found in old GAFF.DEF and new GAFF.DEF Version. But filtering not necessary because of aromatic Ring preclassification to AR2 or AR3
                    if enumToString(mol.atoms.element[prim_neigh]) == "C" && enumToString(mol.atoms.element[sec_neigh]) == "C" 
                        path_stan_stan = string("(", ATD_df.Element_wNeighborCount[prim_neigh],"(", ATD_df.Element_wNeighborCount[sec_neigh],"))")
                        if in(c_tag_list).(path_stan_stan) 
                            c_tag_counter += 1
                        end
                    elseif enumToString(mol.atoms.element[prim_neigh]) == "C" && in(X_dict["XB"]).(enumToString(mol.atoms.element[sec_neigh])) && lastindex(ATD_df.Element_wNeighborCount[sec_neigh]) == "2"
                        path_stan_XB = string("(",ATD_df.Element_wNeighborCount[prim_neigh],"(","XB", ATD_df.Element_wNeighborCount[sec_neigh][lastindex(ATD_df.Element_wNeighborCount[sec_neigh])],"))")
                        if in(c_tag_list).(path_stan_XB)
                            c_tag_counter += 1
                        end
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && lastindex(ATD_df.Element_wNeighborCount[prim_neigh]) == "2" && enumToString(mol.atoms.element[sec_neigh]) == "C"
                        path_XB_stan = string("(","XB", ATD_df.Element_wNeighborCount[prim_neigh][lastindex(ATD_df.Element_wNeighborCount[prim_neigh])],"(", ATD_df.Element_wNeighborCount[sec_neigh],"))")    
                        if in(c_tag_list).(path_XB_stan)
                            c_tag_counter += 1
                        end                            
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && lastindex(ATD_df.Element_wNeighborCount[prim_neigh]) == "2" && in(X_dict["XB"]).(enumToString(mol.atoms.element[sec_neigh])) && lastindex(ATD_df.Element_wNeighborCount[sec_neigh]) == "2"
                        path_XB_XB = string("(","XB", ATD_df.Element_wNeighborCount[prim_neigh][lastindex(ATD_df.Element_wNeighborCount[prim_neigh])],
                                            "(","XB", ATD_df.Element_wNeighborCount[sec_neigh][lastindex(ATD_df.Element_wNeighborCount[sec_neigh])],"))")
                        if in(c_tag_list).(path_XB_XB) 
                            c_tag_counter += 1
                        end                            
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && wgraph_adj[num,prim_neigh] == 1
                            c_tag_counter += 1
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && lastindex(ATD_df.Element_wNeighborCount[prim_neigh]) == "2" && wgraph_adj[num,prim_neigh] == 1
                        c_tag_counter += 1
                    elseif in(X_dict["XD"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && wgraph_adj[prim_neigh,sec_neigh] == 2
                        c_tag_counter += 1
                    end
                end
            end
        end
    end
    if c_tag_counter < 2
        return false
    end
    return true
end


function DEF_e_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "e" tag is added in GAFF.DEF to atomtypes in aromatic non-ring systems with one or more of the e_tag_list properties
    e_tag_list = ["(C2[SB'])", "(C3[SB'])", "(XA1[SB'])", "(XB2[SB'])", "(XD3[SB',db])", "(XD4[SB',db])"]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        if in(ATD_df.BondTypes[prim_neigh]).("NG")
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if in(ATD_df.BondTypes[sec_neigh]).("NG")
                    # build neighbor-"neighbors of neighbor" or bond tuple and check if in c_tag_list
                    if ATD_df.Element_wNeighborCount[prim_neigh] == "C2" && wgraph_adj[num,prim_neigh] == 1
                        return true                            
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && wgraph_adj[num,prim_neigh] == 1
                        return true
                    elseif in(X_dict["XA"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && lastindex(ATD_df.Neighbors[prim_neigh]) == 1
                        return true
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
                        return true
                    elseif in(X_dict["XD"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && wgraph_adj[prim_neigh,sec_neigh] == 2 &&
                        (lastindex(ATD_df.Neighbors[prim_neigh]) == 3 || lastindex(ATD_df.Neighbors[prim_neigh]) == 4)
                        return true
                    end
                end
            end
        end
    end
    return false
end


function DEF_g_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "g" tag is added in GAFF.DEF to atomtypes in aromatic non-ring systems with one or more of the g_tag_list properties
    e_tag_list = ["(C2[SB'])", "(C3[SB'])", "(N1[SB'])", "(XB2[SB'])"]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        if in(ATD_df.BondTypes[prim_neigh]).("NG")
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if in(ATD_df.BondTypes[sec_neigh]).("NG")
                    # build neighbor-"neighbors of neighbor" or bond tuple and check if in c_tag_list
                    if ATD_df.Element_wNeighborCount[prim_neigh] == "C2" && wgraph_adj[num,prim_neigh] == 1
                        return true                            
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && wgraph_adj[num,prim_neigh] == 1
                        return true
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "N1" && wgraph_adj[num,prim_neigh] == 1
                        return true
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
                        return true
                    end
                end
            end
        end
    end
    return false
end


function DEF_x_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "g" tag is added in GAFF.DEF to atomtypes in aromatic non-ring systems with one or more of the g_tag_list properties
    e_tag_list = ["(C3[SB'])", "(XB2[SB'])", "(XD3[sb',db])","(XD4[sb',db])"]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        if in(ATD_df.BondTypes[prim_neigh]).("NG")
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if in(ATD_df.BondTypes[sec_neigh]).("NG")
                    # build neighbor-"neighbors of neighbor" or bond tuple and check if in c_tag_list
                    if ATD_df.Element_wNeighborCount[prim_neigh] == "C2" && wgraph_adj[num,prim_neigh] == 1
                        return true                            
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && wgraph_adj[num,prim_neigh] == 1
                        return true
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && lastindex(ATD_df.Element_wNeighborCount[prim_neigh]) == "2" && wgraph_adj[num,prim_neigh] == 1 && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
                        return true
                    elseif in(X_dict["XD"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && 
                        (lastindex(ATD_df.Neighbors[prim_neigh]) == 3 || lastindex(ATD_df.Neighbors[prim_neigh]) == 4)
                        return true
                    end
                end
            end
        end
    end
    return false
end


function DEF_y_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "g" tag is added in GAFF.DEF to atomtypes in aromatic non-ring systems with one or more of the g_tag_list properties
    e_tag_list = ["(C3[SB'])", "(XB2[SB'])", "(XD3[sb',db])","(XD4[sb',db])"]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        if in(ATD_df.BondTypes[prim_neigh]).("NG")
            for sec_neigh in ATD_df.Secondary_Neighbors[num][i]
                if in(ATD_df.BondTypes[sec_neigh]).("NG")
                    # build neighbor-"neighbors of neighbor" or bond tuple and check if in c_tag_list
                    if ATD_df.Element_wNeighborCount[prim_neigh] == "C2" && wgraph_adj[num,prim_neigh] == 1
                        return true                            
                    elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && wgraph_adj[num,prim_neigh] == 1
                        return true
                    elseif in(X_dict["XB"]).(enumToString(mol.atoms.element[prim_neigh])) && lastindex(ATD_df.Element_wNeighborCount[prim_neigh]) == "2" && wgraph_adj[num,prim_neigh] == 1 && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
                        return true
                    elseif in(X_dict["XD"]).(enumToString(mol.atoms.element[prim_neigh])) && wgraph_adj[num,prim_neigh] == 1 && 
                        (lastindex(ATD_df.Neighbors[prim_neigh]) == 3 || lastindex(ATD_df.Neighbors[prim_neigh]) == 4)
                        return true
                    end
                end
            end
        end
    end
    return false
end


function DEF_h_tag(num::Int64, ATD_df::DataFrame, wgraph_adj::Graphs.SparseMatrix, X_dict::AbstractDict, mol::AbstractMolecule)
    # the "g" tag is added in GAFF.DEF to atomtypes in aromatic non-ring systems with one or more of the h_tag_list properties
    h_tag_list = ["(XX[AR1.AR2.AR3])", "(C3[DB])", "(N2[DB])", "(P2[DB])"]
    for (i,prim_neigh) in enumerate(ATD_df.Neighbors[num])
        # build neighbor-"neighbors of neighbor" or bond tuple and check if in h_tag_list
        if in(X_dict["XX"]).(enumToString(mol.atoms.element[prim_neigh])) && !in(ATD_df.BondTypes[prim_neigh]).("NG") && in(wgraph_adj[prim_neigh,ATD_df.Secondary_Neighbors[num][i]]).(2)
            return true                            
        elseif ATD_df.Element_wNeighborCount[prim_neigh] == "C3" && in(wgraph_adj[prim_neigh,ATD_df.Secondary_Neighbors[num][i]]).(2) && lastindex(ATD_df.Neighbors[prim_neigh]) == 3
            return true
        elseif ATD_df.Element_wNeighborCount[prim_neigh] == "N2" && in(wgraph_adj[prim_neigh,ATD_df.Secondary_Neighbors[num][i]]).(2) && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
            return true
        elseif ATD_df.Element_wNeighborCount[prim_neigh] == "P2" && in(wgraph_adj[prim_neigh,ATD_df.Secondary_Neighbors[num][i]]).(2) && lastindex(ATD_df.Neighbors[prim_neigh]) == 2
            return true
        end
    end
    return false
end

function format_BondTypes!(num::Int64, str_Bonds::AbstractString, ring_class_list::Vector{Vector{String}})
    str_bondtypes_list = Vector{String}()
    non_ring_atom_bool = ("NG" in ring_class_list[num]) 
    for defi = (1:3)
        push!(str_bondtypes_list, enumToString(DefBond(defi)))
    end
    ret_list = Vector{String}()
    for type in str_bondtypes_list
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
        ONSP_present = false
        if true in in(mol.atoms.element[vlist]).([Elements.O,Elements.N,Elements.S,Elements.P])
            ONSP_present = true
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
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) <= 1 && ONSP_present && lastindex(vlist) > 4
            for x in vlist
                push!(ring_class_list[x], "AR2")
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) < 1 && !ONSP_present && lastindex(vlist) > 4
            # check if Ring vlist has intersections with other rings in molecule and if these are aromatic
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
        elseif (pi_elec / lastindex(vlist)) < 1/2
            for x in vlist
                if !in(ring_class_list[x]).("AR5") && !in(ring_class_list[x]).("AR1")
                    push!(ring_class_list[x], "AR5")    
                end
            end
        end
    end
    return ring_class_list
end


function count_EWG(num::Int64, ATD_df::DataFrame)
    # Electron withdrawal Atoms according to antechamber document are only: N, O, F, Cl, and Br
    # To Do: Test differences for Atom Typing, see below typical know EWG
    strong_pullers = ["Cl1", "F1", "Br1", "I1", "O1", "S1"]
    possible_pullers = ["C2", "C3", "C4", "S3", "N3", "P3", "P4", "O2", "S2"]
    elec_pullers_num = 0
    if true in in(ATD_df.Element_wNeighborCount[ATD_df.Neighbors[num]]).(possible_pullers)
        for i = (1:lastindex(ATD_df.Neighbors[num]))
            if true in in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]])
                # all neighbors of neighbor that are in strong_pullers are an EWG
                elec_pullers_num += countmap(in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]]))[true]
            end
        end
    end
    if elec_pullers_num == 0
        elec_pullers_num = -1
    end
    # # println("EWG: ", elec_pullers_num)
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
        # println("no cycle detected")
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