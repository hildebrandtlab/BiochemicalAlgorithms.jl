using BiochemicalAlgorithms

export load_all, run_atomtyping, string_as_variable_name, conti_return, testing_loop_conti

function load_all()
    file_location = "test/data/gaff_paper_examples/a/"
    mol_df = DataFrame([Vector{String}(undef, lastindex(readdir(file_location))), Vector{AbstractMolecule}(undef, lastindex(readdir(file_location)))], 
                        ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location))
        mol_df.molname[num] = string("mol", i[1:2])
        mol_df.abstract_mol[num] = load_pubchem_json(string(file_location, i))
    end
    return mol_df
end

function run_atomtyping()
    mol_df = load_all()
    df = select_atomtyping()
    exit_dict = Dict{Symbol, DataFrame}
    for num = (1:nrow(mol_df))
        println(string(mol_df.molname[num]))
        atomtypes_list = get_atomtype(mol_df.abstract_mol[num], df)
        exit_dict = merge(exit_dict, Dict(Symbol(string(mol_df.molname[num],"atomtypes")) => atomtypes_list))
    end
    return exit_dict
end

function string_as_variable_name(str::AbstractString, var::Any)
    str = Symbol(str)
    return @eval (($str) = ($var))
end
