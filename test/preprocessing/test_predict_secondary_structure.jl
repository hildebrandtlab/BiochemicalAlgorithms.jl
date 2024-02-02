

@testitem "predict secondary structure for BPTI" begin

    default_fdb = FragmentDB()
    fdb = FragmentDB(ball_data_path("fragments/Fragments.db.json"))
    sys = load_pdb(ball_data_path("../test/data/PDBFile_test2.pdb"))
   
    normalize_names!(sys, fdb)
    reconstruct_fragments!(sys, fdb)
    build_bonds!(sys, fdb)
    predict_hbonds!(sys, :KABSCH_SANDER)
    summary = predict_secondary_structure!(sys)
    
    #@test natoms(sys) == 892 ## vorher 892 TODO checken ob das an Cysteinen liegt

    # the reference was generated by the C++ BALL version
    @test String(summary) == "--GGGG-----------EEEEEEETTTTEEEEEEE---------B--HHHHHHHH---"
    @test length(secondary_structures(sys)) == 9
    
    # the following was additionally 'compared to' in the C++ version 
	#String PDB_summary = "CCGGGGSCCCCCSCCCCEEEEEEETTTTEEEEEEECSSSCCSSCBSSHHHHHHHHSCC";
    


    function predict_secondary_structure(sys)
        default_fdb = FragmentDB()
        fdb = FragmentDB(ball_data_path("fragments/Fragments.db.json"))

        filter_fn = a -> is_amino_acid(parent_fragment(atom_by_idx(sys, a.idx)))

        filtered_sys = filter_atoms(
            filter_fn,
            sys
        )
        
        normalize_names!(filtered_sys, fdb)
        reconstruct_fragments!(filtered_sys, fdb)
        build_bonds!(filtered_sys, fdb)
        predict_hbonds!(filtered_sys, :KABSCH_SANDER)
        summary = predict_secondary_structure!(filtered_sys)

        return summary
    end

    ### 5pti.pdb
    sys_5pti = load_pdb(ball_data_path("../test/data/5PTI.pdb"))
    @test String(predict_secondary_structure(sys_5pti)) == 
          "--GGGG-----------EEEEEEETTTTEEEEEEE---------B--HHHHHHHH---"
    
    #### 3fgu.pdb
    sys_3fgu = load_pdb(ball_data_path("../test/data/3fgu.pdb"))
    @test String(predict_secondary_structure(sys_3fgu)) ==
          "-HHHHHHHHHHHHHHHGGG---HHHHHHHHHHHHHHHHHHH-TTTTT------EE---B------EEEEEEEE---EEEEEEEEEETTTEEEE-EEEEEE--HHHHH-BHHHHHHHHHHHHHHHHHTTT----EEEEEE--EEEE-BTTEEEE-B--TT---BT-BT-BHHHHHHHHHHHH-----EE-EEE-HHHHHHHHTTTT-TTB-EEEEE---EEEEEEEEGGG-TT------EEEEE--GGGTTTT-TTGGG--HHHHHHHHT-TTTT--TTGGG--HHHHHHHHHHHHHHHHHTT--GGG---TTTT-TT---HHHHHHHHT-----HHHHHHHHHTT----HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTT---EEEEEEEE-HHHHH---HHHHHHHHHHHH-TTEEEEEEE---HHHHHHHHHHHHH--"

    ### 2ptc.pdb
    #sys_2ptc = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
    #println("2ptc:", String(predict_secondary_structure(sys_2ptc)))
    #println("    :", "-BT-EE--TT--TTEEEEE---B-EEEEE-BTTEEEE-GGG-----EEEE----TT------EEEEE-EEEE-TT--TTT-TT--EEEEE-------------EE------TT-EEEEEE--B-------B----EEEEEEE--HHHHHHHTTT---TTEEEEE-TT---B--TT-TT-EEEETTEE-EEE-B------TT--EEEEEGGGGHHHHHHHHHT-")

    ### 6nxl.pdb
    #sys_6nxl = load_pdb(ball_data_path("../test/data/6nxl.pdb"))
    #println("6nxl:", String(predict_secondary_structure(sys_6nxl)))
    #println("    :", "-------------------TT-BHHHHHHHHHHHH---GGGEEEEETTEE--TT-BTGGGT--TT---EEEE---")

    #println(fragments_df(sys).name)

end