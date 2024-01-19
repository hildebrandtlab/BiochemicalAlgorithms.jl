export predict_secondary_structure!

@enum _DSSP_bridge_pattern PARALLEL ANTIPARALLEL DEFAULT

DSSP_DEBUG = true #false

function _insert_turn!(target_n_turn::String, turn_char::Char, position::Int)
	# first sign
	if target_n_turn[position] ∈ ['<', 'X']
		target_n_turn[position] = 'X';
	else
		target_n_turn[position] = '>';
	end
	#  positions in between
	for j in 1:turn-1
		if target_n_turn[position+j] ∈ ['-']
			target_n_turn[position+j] = turn_char
		end
	end
	# last position 
	if target_n_turn[position+turn] ∈ [ '>' ,'X']
		target_n_turn[position+turn] = 'X';
	else
		target_n_turn[position+turn] = '<';
	end	
end

function _test_pattern_at_two_positions(target_string, pos_i::Int, test_char::Char)
    #function test_string(const String& s, Size offset, Size offset_2)
        #testString2_(s, offset) && testString2_(s, offset + offset_2)
    return  _test_pattern_at_position(target_string, pos_i)  &&
            _test_for_char_at_position(target_string, pos_i, test_char)
end
    
function _test_pattern_at_position(target_string, pos_i::Int)
        #if (s.size() < offset + 2) return false;
        if length(target_string) <= pos_i + 1 # TODO +1??
            return false
        end

        return target_string[pos_i:pos_i+1] ∈ [">>", "XX", ">X", "X>"]
end
    
function _test_for_char_at_position(target_string, pos_i::Int, test_char::Char)
    #if (s.size() < offset + 2) return false;
    if length(target_string) <= pos_i + 1 # TODO +1??
        return false
    end

    if target_string[pos_i:pos_i+1] ∈ ["><", "X<"]
        return true
    end

    return (( target_string[pos_i+1] == test_char )     
        &&  (target_string[pos_i] ∈ [ '>', 'X'])) 
end
    
function _prepare_hbonds(ac::AbstractAtomContainer{T}) where {T}
    raw_hbonds = backbone_hydrogen_bonds(ac)

    #if DSSP_DEBUG
    #    println("DSSP_DEBUG raw_hbonds: ", raw_hbonds)
    #end

    # we cannot rely on residue numbers or indices, these might not be continuous
    # e.g., after removing water
    # thus, we map fragment indices to the position of the fragment in the fragments-Vector
    # of the given AtomContainer
    idx_map = Dict(f.idx => i for (i, f) in enumerate(fragments(ac)))

    # the hbonds-vector contains the bonds, i.e., connects atoms; we need to map this
    # to pairs of fragment indices (needed for DSSP)
    # these are stored such that the acceptor is mapped to the donor

    function convert_bond(b::Bond)
        a1 = atom_by_idx(ac, b.a1)
        a2 = atom_by_idx(ac, b.a2)

        if a1.name == "N"
            (idx_map[parent_fragment(a2).idx], idx_map[parent_fragment(a1).idx])
        else
            (idx_map[parent_fragment(a1).idx], idx_map[parent_fragment(a2).idx])
        end
    end

    bonded_fragments = convert_bond.(raw_hbonds)

    fragment_list = [findall(t -> t[1] == i, bonded_fragments) for i in 1:nfragments(ac)]

    result = [map(t->t[2], d) for d in getindex.(Ref(bonded_fragments), fragment_list)]
    if DSSP_DEBUG
        println("DSSP_DEBUG converted hbonds: ", result)
        for (i, f) in enumerate(result)
            println(i,":", f)  
        end  
    end

    result
end


"""
    $(TYPEDSIGNATURES)

    Predict secondary structure from 3D structure for a given AtomContainer.

    The implementation follows the DSSP algorithm described in
        "Kabsch W & Sander C (1983). 
        Dictionary of protein secondary structure: 
        pattern recognition of hydrogen-bonded and geometrical features. 
        Biopolymers, 22, 2577-2637."

    TODO: When applied to a protein, it removes the instances of SecondaryStructure
    from the protein, predicts the secondary structure elements based
    (mostly) on H-bond patterns and backbone torsions and reinserts the
    appropriate secondary structure elements at the predicted positions.

    A previous H-Bond prediction for the backbone is needed, as implemented in  
    the KABSCH_SANDER variant of the predict_hbonds! - method
    


    Example:

    ```julia
    fdb = FragmentDB()
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)
    predict_hbonds!(sys, :KABSCH_SANDER)
    predict_secondary_structure!(sys)
    ```
"""
function predict_secondary_structure!(ac::AbstractAtomContainer{T}) where {T<:Real}
    # the following should happen outside
    # fdb = FragmentDB()
    # normalize_names!(ac, fdb)
    # reconstruct_fragments!(ac, fdb)
    # build_bonds!(ac, fdb)
    # predict_hbonds!(sys, :KABSCH_SANDER)

    # get the backbone h_bond pattern according to KABSCH_SANDER
    hbonds = _prepare_hbonds(ac)

    #   initialize the summary strings	
    size = length(hbonds)  
    
    sheet     = Vector{Char}(  repeat('-', size))
    fiveturn  = Vector{Char}(  repeat('-', size))
    fourturn  = Vector{Char}(  repeat('-', size))
    threeturn = Vector{Char}(  repeat('-', size))
    summary   = Vector{Char}(  repeat('-', size))


    # note that hbonds.size() is not the number of hbonds but
    # rather the number of residues in the system.

    # first search turns
    for i in eachindex(hbonds) 
        

        # over all HBondpartners
        for k in 1:length(hbonds[i])
       
            # ---------- 4 turns ----------
            #    >444<
            if hbonds[i][k] == (i+4)

                # first position
                #if  (fourturn[i]  == '<') || (fourturn[i] == 'X')
                if  fourturn[i] ∈ ['<', 'X']
                    fourturn[i]   = 'X';
                else
                    fourturn[i]   = '>'
                end
                # position 2,3,4
                for j in 1:3 ## (int j=1; j<4; j++)
                    if fourturn[i+j] == '-'
                        fourturn[i+j] = '4'
                    end
                end

                # last position has always to be <
                fourturn[i+4] = '<'
            # ---------- 3 turns -------------
            elseif hbonds[i][k] == (i + 3)

                #println("i:", i, " a threeturn is detected:", k)

                #if threeturn[i] == '<' || threeturn[i] == 'X'
                if threeturn[i] ∈ ['<', 'X']
                    threeturn[i] = 'X'
                else
                    threeturn[i] = '>'
                end  
                # the middle part
                for j in 1:2
                    if threeturn[i+j] == '-'
                        threeturn[i+j] = '3'
                    end
                end
                # last position has always to be <
                threeturn[i+3] = '<'

            # ---------- 5 turns ----------
            elseif hbonds[i][k] == (i + 5)
                #if fiveturn[i] == '<' || fiveturn[i] == 'X'
                if fiveturn[i] ∈ ['<', 'X']
                    fiveturn[i] = 'X'
                else
                    fiveturn[i] = '>'
                end
                # the middle part
                for j in 1:4
                    if fiveturn[i+j] == '-'
                        fiveturn[i+j] = '5'
                    end
                end
                #last position has always to be <
                fiveturn[i+5] = '<'
            end

        end
        
    end

    # ****************************************
    # * in the next step, we search bridges
    # ****************************************
    #initialize posbridges 
    posbridges = Vector{Vector{Int}}(repeat([[]], size)) 
    
    # iterate over all 'residues' 
    for current_res in 1:length(hbonds) 
        # nothing found => nothing done
        if length(hbonds[current_res]) == 0 
            continue
        end

        # over all HBondpartners
        for current_bond in 1:length(hbonds[current_res])  	
          
            partner = hbonds[current_res][current_bond]

            # nothing found => nothing done
            if partner == 0
                continue
            end

            # actually we have HBP(i,k) with  i=current_res, k=partner
            # do we have HBP(k, i+2) ?  => parallel bridge(i+1,k) 			
            for s in 1:length(hbonds[partner])
              
                if hbonds[partner][s] == current_res + 2
                  
                    push!(posbridges[current_res+1], partner)
                    
                    push!(posbridges[partner], current_res+1)
                 
                    sheet[current_res+1] = '+'
                 
                    sheet[partner] = '+'
                end
               
            end # forend parallelbridge
  
            # currently we have HBP(i,k) with i=current_res and k=partner
            # do we have HBP(k, i) ? => antiparallel bridge(i,k)
            for s in 1:length(hbonds[partner])
                if hbonds[partner][s] == current_res
                    # we aren't allowed to overwrite antiparallel bridges found before!
                    # remember: we have two equal cases: 
                    # 	 	first : HBP(k,i) && HBP(i,k)   and
                    # 	 	second: HBP(i,k) && HBP(k,i)   						
                        
                    # insert
                    # NOTE: there might be  more than two bridges for this residue
                    # but this should never happen ;-)
                    push!(posbridges[current_res], partner)
                    push!(posbridges[partner], current_res)
                    sheet[current_res] = '/'
                    sheet[partner] = '/'

                end	
            end #forend antiparallelbridge

          

            # currently we have HBP(i,k) with i=current_res and k=partner
            # do we have HBP(k-2, i+2) ? => antiparallel bridge(i+1,k-1)
            ###if (partner-2 >= 0) && (current_res+2 <= size)   size = length(hbonds)
            if (partner-2 >= 1) && (current_res+2 <= length(hbonds))
                for s in 1:length(hbonds[partner-2])
                    if (hbonds[partner-2][s] == current_res + 2) && (current_res + 1 != partner - 1)
                        # insert
                        # NOTE: there might be more than two bridges for this residue
                        # but this should never happen ;-)
                
                        push!(posbridges[current_res+1], partner-1)
                        push!(posbridges[partner-1], current_res+1)
                        sheet[current_res+1] = '/'
                        sheet[partner-1] = '/'
                    end
                end

            end # if (partner-2>=...)
        end
     
    end

   
    #
    #  now we search for ladders!
    #
    # ladder: set of one or more consecutive bridges of identical type
    # or 
    # bulge-linked ladders: two ladders or bridges of the same type 
    #                          connected by at most one extra residue 
    #                          on one strand and at most four extra 
    #                          residues on the other strand
    
    if DSSP_DEBUG
        println("DSSP_DEBUG: now we search for ladders")
    end

    parallel     = 'a' - 1
    antiparallel = 'A' - 1
    letter       = '-'

    last_pattern    = DEFAULT
    found_a_pattern = false
    no_residue      = -5
    last_parallel_res     = no_residue
    last_antiparallel_res = no_residue
    last_residue          = no_residue

    #over all residues
    for residue in 1:length(hbonds) 
        if residue > length(posbridges)
            break
        end
        
        # do we have a bridge?		
        if sheet[residue] != '-'
            # parallel bridge	
            if sheet[residue] == '+' || (sheet[residue] > ('a'-1) && sheet[residue] < ('z'+1))
                last_residue = last_parallel_res
                last_pattern = PARALLEL
            else #  antiparallel bridge
                last_residue = last_antiparallel_res
                last_pattern = ANTIPARALLEL
            end
            found_a_pattern = false
            has_bonds_to_the_left = false

            if last_residue != no_residue
                for last_part_i in 1:length(posbridges[last_residue])
                    for curr_part_i in 1:length(posbridges[residue])
                        last_part = posbridges[last_residue][last_part_i]
                        curr_part = posbridges[residue][curr_part_i]

                        # have we already seen this residue by a bridge before?
                        if curr_part < convert(Int, residue)
                            has_bonds_to_the_left = true
                            letter = sheet[curr_part]
                        end

                        # do we have a continuation of a ladder?
                        #   allowed are bulge-linked ladders 
                        #   which consist of two ladders or bridges of the same type 
                        #   connected by at most one extra residue on one strand and 
                        #   at most four extra residues on the other strand!   
                        r_diff = residue - last_residue;
                        p_diff = abs(curr_part - last_part);

                        if ((r_diff == 1 && p_diff == 1) ||
                            (r_diff  < 2 && p_diff  < 5) ||
                            (r_diff  < 5 && p_diff  < 2)
                        )
                            found_a_pattern = true
                        end
                    end
                end
            end


            # Do the naming
            # NOTE: there is no problem if we have seen the 
            # residue before and set the letter 
            #  AND found a pattern and overwrite the letter perhaps! 
            #  the sheet - Loop(see below) will give them a unique letter 
            if (last_residue != no_residue && found_a_pattern)
                letter = sheet[last_residue];	
            else 
                if !has_bonds_to_the_left
                    if last_pattern == PARALLEL
                        parallel += 1;
                        if parallel == 'z'+1 
                            parallel = 'a'
                        end
                        letter = parallel;	
                    else 
                        antiparallel+=1;
                        if (antiparallel == 'Z'+1) 
                            antiparallel='A';
                        end
                        letter = antiparallel;	
                    end				
                end
            end
            if last_pattern == PARALLEL
                last_parallel_res = residue;
            else
                last_antiparallel_res = residue;
            end
            
            # name all residues belonging to this bridge
            # NOTE: in most cases there will be just one entry
            sheet[residue] = letter;	

            ##### NOTE:   SKIPPED by amoll for unknown reason
            #for  curr_part in 1:length(posbridges) 
            #    if posbridges[residue][curr_part] > convert(Int, residue)
            #    #if we have a bridge to the right, we name the partner
            #        #	sheet[posbridges[residue][curr_part]] = letter;
            #    end
            #    #		sheet[posbridges[residue][curr_part]] = letter;
            #end
            #### END SKIPPED
        end
    end # end of search ladders

    if DSSP_DEBUG
        println("DSSP_DEBUG: naming (looking for ) sheet_s")
    end 

    #
    # now we are naming (looking for ) sheet_s
    #
    # sheet: set of one or more ladders connected by shared residues
    for residue in 1:length(hbonds) 
        if residue > length(posbridges)
            break
        end

        # do we have a bridge?		 	
        if (sheet[residue] != '-')
            letter = sheet[residue];
            for curr_part_i in 1:length(posbridges[residue])
            
                curr_part = posbridges[residue][curr_part_i];			
                if sheet[curr_part] != letter
                    ##changeAllXToY_(sheet[curr_part], letter, sheet); 
                    sheet = replace(sheet, sheet[curr_part] => letter)
                end				
            end				
        end
    end # end of sheet_s

    if DSSP_DEBUG
        println("DSSP_DEBUG: now we repair the irregularities")
    end

    # ***********************************
    # * now we repair the irregularities
    # ***********************************
    #  two overlapping minimal helices(consecutive n-turns)
    #  offset by two or three residues are joined into one 
    #  helix
    #  therefore we just insert the missing n-turns

    turn_length = 5;
    #turn = 3;

    # we start with the 3-turn
    #  -------------- 3-turn -----------------------
    # we have to check the length : offset by two or three 
    # residues and the length of the turn itself (here 5)
    for i in 1:length(hbonds)
  
        # case offset two residues 			
        if (i + 2 + turn_length <= length(hbonds)	&&
                _test_pattern_at_two_positions(threeturn, i, '3') &&  
                fourturn[i + 2] ∉ ['>', 'X'] 			
            )
            _insert_turn!(threeturn, '3', i+2) ## TODO offby1 ?
        end

        # case offset three residues	
        if (i + 3 + turn_length <= length(hbonds)     &&  
            _test_pattern_at_two_positions(threeturn, i, '4') &&
            fourturn[i+2] ∉ [ '>', 'X'] &&
            fourturn[i+3] ∉ [ '>', 'X'] 
            )
            _insert_turn!(threeturn, '3', i+2) ## TODO offby1
            _insert_turn!(threeturn, '3', i+3) ## TODO offby1
        end		
    end

    #  -------------- 4-turn -----------------------
    # we have to check the length : offset by two or three 
    # residues and the length of the turn itself (here 6) 
    turn_length = 6;
    #turn = 4;	
    for i in 1:length(hbonds)
        # case offset two residues
        if (i + 2 + turn_length <= length(hbonds)	&&
            _test_pattern_at_two_positions(fourturn, i, '3') 	&&
            fourturn[i + 2]  ∉ ['>', 'X']
        )
            _insert_turn!(fourturn, '4', i+2) ## TODO offby1 ?
        end		

        # case offset three turns	
        if (i + 3 + turn_length <= length(hbonds) 	&&  
            _test_pattern_at_two_positions(fourturn, i, '4') 	&&
            fourturn[i + 2]  ∉ ['>', 'X'] &&
            fourturn[i + 3]  ∉ ['>', 'X']		
        )
            _insert_turn!(fourturn, '4', i+2) ## TODO offby1 ?
            _insert_turn!(fourturn, '4', i+3) ## TODO offby1 ?	
        end		
    end

    # -------------- 5-turn -----------------------
    # we have to proof the length : offset by two or three 
    # residues and the lenght of the turn itself (here 7) 

    turn_length = 7;
    #turn = 5;			

    for i in 1:length(hbonds)
        # case  offset two residues
        if (i + 2 + turn_length  <= length(hbonds)	&& 
            _test_pattern_at_two_positions(fiveturn, i, '3') 	&&    
            fourturn[i + 2]  ∉ ['>', 'X']	
        )
            _insert_turn!(fiveturn, '5', i+2) ## TODO offby1 ?
        end		

        # case offset three turns	
        if (i + 3 + turn_length  <= length(hbonds)	&&
            _test_pattern_at_two_positions(fiveturn, i, '4') &&
            fourturn[i + 2]  ∉ ['>', 'X']  &&
            fourturn[i + 3]  ∉ ['>', 'X']
        )
            _insert_turn!(fiveturn, '5', i+2) ## TODO offby1 ?
            _insert_turn!(fiveturn, '5', i+3) ## TODO offby1 ?
        end		
    end
    
    # *****************************************************
    # *
    # *     now we construct the summary string  
    # *
    # * structural overlaps are eleminated by considering
    # * hierarchy H > B > E > G > I > T
    # * 		H means 4 Helices, 
    # *         B means single bridges, 
    # * 		E means extended ladder residues, 
    # *         G means 3 Helices
    # *		    I means 5 Helices and 
    # *         T means single 3-, 4-, or 5-turns 
    # * we start with writing 5 Helices and overwrite graduately 
    # * the summary_ string with 3 Helices, extended bridges, 
    # * single bridges and 4 Helices helices	
    # *
    # ****************************************************
    

    # --------------------- 5 helices ---------------------- 
    for i in 1:length(hbonds) 
        # we initialize the summary_ string with '-'			
        if (fiveturn[i] == '-') 
            summary[i] = '-';
        #else if (testString2_(fiveturn, i))
        elseif _test_pattern_at_position(fiveturn, i)
            if i+5 < length(hbonds) ## TODO i+5 <= length(hbonds)?
                summary[i+1] = 'I';
                summary[i+2] = 'I';
                summary[i+3] = 'I';
                summary[i+4] = 'I';
                summary[i+5] = 'I';
            end	
        else   #  do we have a helix reduced to less than minimal size?
            #String ss = fiveturn_.getSubstring(i).toString();
            #if (testString3_(fiveturn_, i, '5'))
            if _test_for_char_at_position(fiveturn, i, '5') 
                #for (int j=1; (j<5) && ((i+j)<summary_.size()) && ((i+j)<fiveturn_.size()) ;j++)
                for j in 1:5
                    if ((i+j) > length(summary)) ||  ((i+j) > length(fiveturn))
                        break
                    end
                    if (fiveturn[i+j] != 'I')
                        summary[i+j] = 'T';
                    end				
                end
            end	
        end
    end	

    

    # -------------------3 helices ------------------------
    for i in 1:length(hbonds)
        #if (testString2_(threeturn_, i))
        if _test_pattern_at_position(threeturn, i)
            if ( (i+3) < length(hbonds)) # TODO?? i+3 <= length(hbonds)?
                summary[i+1]= 'G';
                summary[i+2]= 'G';
                summary[i+3]= 'G';
            end	
        # do we have a helix reduced to less than minimal size?
        # we have to consider, that we do not overwrite 
        #else if(testString3_(threeturn_, i, '3'))
        elseif _test_for_char_at_position(threeturn, i, '3') 
            if i+3 < length(summary) # TODO?? i+3 <= length(hbonds)?
                for j in 1:3 #(Size j=1; j<3;j++)
                    if summary[i+j] ∉ ['G',  'I']
                        summary[i+j] = 'T';	
                    end
                end
            end	
        end	
    end


    # ---------------- Extended Bridges and Single Bridges --------------
    # according to the paper:
    # 		single bridges are ladders of length one -> B, 
    # 		all other ladder residues -> E
    # we assume that there is a mistake in the paper: E has a higher priority than B
    # first we generate the sheet_-line and than summarize it in the summary_ line

    #for(Size i=0; i< (sheet_.size()); i++)
    sheet_test_index = 1
    while sheet_test_index < length(sheet)
        if sheet[sheet_test_index] != '-'
            letter = sheet[sheet_test_index]
            
            # check the length of that sheet
            j::Int = 0 # start a sheet	
            #for(j=0; ((i+j)<sheet_.size()) && (sheet_[i+j]==letter); j++)
            while ((sheet_test_index+j) <= length(sheet) && (sheet[sheet_test_index+j]==letter))
                j += 1
            end
            # case length = 1
            if ((j==1) && (summary[sheet_test_index] != 'E')) #single bridge
                summary[sheet_test_index] = 'B'
            #elseif j==0 ## should not happen!!!
            #    #println("uups")
            else #extended bridge
                for n in 0:j-1	
                #for(int n=0; n<j; n++)
                    summary[sheet_test_index + n] = 'E';
                end
                # hop forward
                sheet_test_index = sheet_test_index + j - 1;	
            end
        end		
        sheet_test_index += 1								
    end


    #				
    #  ---------------- 4 helices ---------------------
    # 
    #for(Size i= 0; i<size; i++)
    for i in 1:length(hbonds)
        #if (testString2_(fourturn_, i))
        if _test_pattern_at_position(fourturn, i)
            if (i+4 < length(hbonds))
                summary_[i+1] = 'H';
                summary_[i+2] = 'H';
                summary_[i+3] = 'H';
                summary_[i+4] = 'H';
            end	
        # or do we have a helix reduced to less than minimal size?
        # we have to ensure, that we do not overwrite 
        elseif (i < length(hbonds) ) && (fourturn[i:i+1] ∈ [['>','4'], ['X','4']]  ) 
            # single 3-, 4- or 5- helix
            if i+4 < length(hbonds)
                #for(Size j=1; j<4; j++)
                for j in 1:3
                    # do not overwrite!
                    if summary[i+j] ∉ ['G','H','I','E','B']
                        summary[i+j]= 'T';	
                    end
                end
            end #if	
        end #if	
    end

    # finally we need to check the summary string again 
    # in order to identify and correct 'single/shortend' G or Is, 
    # (i.e. remains generated by partial overwrite a pattern GGG or IIIII by HHHH )

    #for(Size i=0; i<( summary_.size()); i++)
    for i in 1:length(summary)
                
        if  (i+2) <= length(summary)   
            if ((  (summary[i]   != 'G')
                && (summary[i+1] == 'G')
                && (summary[i+2] != 'G'))
            || (   (summary[i]   != 'I')
                && (summary[i+1] == 'I')
                && (summary[i+2] != 'I'))
            )    
                # [!G,!H][G,H][!G,!H] => [!G,!H] T [!G,!H]
                summary[i+1] = 'T';
            end
        end

        if  (i+3) <= length(summary) 
            if ((  (summary[i]   != 'G')
                && (summary[i+1] == 'G')
                && (summary[i+2] == 'G')
                && (summary[i+3] != 'G')) 
            || (   (summary[i]   != 'I')
                && (summary[i+1] == 'I')
                && (summary[i+2] == 'I')
                && (summary[i+3] != 'I'))
            )
                summary[i+1]='T';
                summary[i+2]='T';
            end
        end	
        
        if (     ((i+4) <= length(summary))
            &&  (summary[i]  !='I')
            &&  (summary[i+1]=='I')
            &&  (summary[i+2]=='I')
            &&  (summary[i+3]=='I')
            &&  (summary[i+4]!='I')
        )				
            summary[i+1]='T';
            summary[i+2]='T';
            summary[i+3]='T';
        end	
    end

    if DSSP_DEBUG
        println("DSSP DEBUG resulting summary:   ", String(summary))
        println("DSSP DEBUG resulting threeturn: ",String(threeturn))
        println("DSSP DEBUG resulting fourturn:  ", String(fourturn))
        println("DSSP DEBUG resulting fiveturn:  ", String(fiveturn))
        println("DSSP DEBUG resulting sheet:     ", String(sheet))
        println("DSSP DEBUG resulting summary:   ", String(summary))
    end

    # DONE DSSP!

    return summary
end