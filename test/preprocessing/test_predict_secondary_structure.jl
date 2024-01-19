

@testitem "predict secondary structure for BPTI" begin

    function get_unified_secondary_structure_type(type)
        result = 'L'
        if type == 'H'  # ss->setType(SecondaryStructure::HELIX);
						# ss->setProperty("HELIX_TYPE", "ALPHA");
			result = 'H'
		elseif type == 'G'  # ss->setType(SecondaryStructure::HELIX);
						    # ss->setProperty("HELIX_TYPE", "THREE_TEN");
			result = 'G'
        elseif type == 'I'  # ss->setType(SecondaryStructure::HELIX);
						    # ss->setProperty("HELIX_TYPE", "PI");
			result = 'I'
	    elseif type == 'E'  # ss->setType(SecondaryStructure::STRAND);
			result = 'E';				 
        else  #ss->setType(SecondaryStructure::COIL);
			result = 'L'
        end
        #println("DEBUG TEST get_unified_secondary_structure_type(type)" , type, "--->", result)
        return result
	end

    function get_secondary_structures(summary)
         
        #println("DEBUG DSSP test get_secondary_structures: given ", summary)

        # 
        last_ss_type = 'X'
        #println(last_ss_type)
         
        # each tuple contains
        #  first:  the ss type
        #  second: the index of the starting residue
        secondary_structures = Vector{Tuple{Char,Int64}}() 
        
        #println("laenge summary:", length(summary))
        
        for i in 1:length(summary)
            
            #println("DEBUG DSSP test get_secondary_structures: i ", i, "last_type =",  last_ss_type , ",", summary[i])
            
            # a switch in the SS type
            # depending on the last type of secondary structure we have seen,
		    # we need to react differently to merge them sensibly
            if last_ss_type != summary[i]
                if last_ss_type == 'L' 
                    # we are in a loop
                    # note that we identify 'real' loops, isolated bridges and turns (-,B,T)
                    # and map them all to loops. Thus we need to determine here if the current
                    # residue also maps to a loop (we already know that the last residue was one)
                    if summary[i] ∈ ['-', 'B','T']
                        continue; # nothing to see here... please walk on...
                    else # the current residue is not of type loop => build a new SecondaryStructure 
                        last_ss_type = get_unified_secondary_structure_type(summary[i])
                        ss = (last_ss_type, i)
                        push!(secondary_structures, ss)

                        println("                 new ss element l case:", summary[i], last_ss_type, i)
                    end
                else	# in all other cases, setSecondaryStructure does the hard work
                    last_ss_type = get_unified_secondary_structure_type(summary[i])
                    ss = (last_ss_type, i)
                    push!(secondary_structures, ss)

                    println("                 new ss element else", summary[i], last_ss_type, i)
                end     
            end 
            #resnum += 1
        end
        #println("TEST RESULT:")
        #println(summary)
        #println(secondary_structures)
        
        secondary_structures
    end

    default_fdb = FragmentDB()
    fdb = FragmentDB(ball_data_path("fragments/Fragments.db.json"))
    sys = load_pdb(ball_data_path("../test/data/PDBFile_test2.pdb"))
   
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)
    predict_hbonds!(sys, :KABSCH_SANDER)
    summary = predict_secondary_structure!(sys)
 
    println("BPTI test result: ", String(summary))

    println("soll            : --GGGG-----------EEEEEEETTTTEEEEEEE---------B--HHHHHHHH---")


   
    #summary = ['-', '-', '-', '-', '-', 'H', 'H', 'H', 'H', '-', 
    #           'L', '-', 'B', '-', 'T', '-', '-', '-', '-', '-', 
    #           '-', 'T', 'T', 'T', '-', '-', '-', '-', '-', '-', 
    #           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 
    #           'E', 'E', 'E', 'E', 'E', '-', '-', '-', '-', '-', 
    #           '-', '-', '-', '-', '-', '-', '-', '-']
    
    @test length(get_secondary_structures(summary)) == 7
	#String PDB_summary = "CCGGGGSCCCCCSCCCCEEEEEEETTTTEEEEEEECSSSCCSSCBSSHHHHHHHHSCC";
    
    #@test natoms(sys) == 898 ## vorher 892 TODO checken ob das an Cysteinen liegt
    
    println(fragments_df(sys).name)

end