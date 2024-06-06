fileformat_pdb = SUITE["FileFormats"]["PDB"] 

fileformat_pdb["Reading"] = @benchmarkable load_pdb(ball_data_path("../benchmark/data/AmberFF_bench.pdb")) 


fileformat_json = SUITE["FileFormats"]["Json"] 
fileformat_json["Reading"] = @benchmarkable load_pubchem_json(ball_data_path("../benchmark/data/aspirin_pug_bonds.json")) 


fileformat_json = SUITE["FileFormats"]["SDfile"] 
fileformat_json["Reading"] = @benchmarkable load_sdfile(ball_data_path("../test/data/sdfile_test_1.sdf"))

#s = ball_data_path("../test/data/sdfile_test_1.sdf")
#fileformat_json["Writing"] = @benchmarkable write_sdfile(s, ball_data_path("../benchmark/data/sdf_out.sdf"))
