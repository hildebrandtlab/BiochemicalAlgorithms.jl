using BiochemicalAlgorithms

export load_all

function load_all()
    file_location = "C:/Users/samhu/source/repos/BiochemicalAlgorithms.jl/test/data/gaff_paper_examples/a/"
    mol_list = DataFrame([Vector{String}(undef, lastindex(readdir(file_location))), Vector{AbstractMolecule}(undef, lastindex(readdir(file_location)))], ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location))
        mol_list.molname = string("mol", i[1:5])
        mol_list.abstract_mol[num] = BiochemicalAlgorithms.load_pubchem_json(string(file_location, i))
    end
    println(mol_list)
end