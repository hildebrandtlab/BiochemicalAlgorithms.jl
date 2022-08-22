using BiochemicalAlgorithms

export load_all, run_atomtyping, string_as_variable_name

function load_all()
    file_location = "C:/Users/samhu/source/repos/BiochemicalAlgorithms.jl/test/data/gaff_paper_examples/a/"
    mol_df = DataFrame([Vector{String}(undef, lastindex(readdir(file_location))), Vector{AbstractMolecule}(undef, lastindex(readdir(file_location)))], 
                        ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location))
        mol_df.molname[num] = string("mol", i[1:3])
        mol_df.abstract_mol[num] = BiochemicalAlgorithms.load_pubchem_json(string(file_location, i))
    end
    return mol_df
end

function run_atomtyping()
    mol_df = load_all()
    df = select_atomtyping()
    for num = (1:nrow(mol_df))
        # string_as_variable_name(mol_df.molname[num], mol_df.abstract_mol[num])
        curr_mol = mol_df.abstract_mol[num]
        curr_mol_name = Symbol(mol_df.molname[num])
        atomtypes_list = get_atomtype(mol_df.abstract_mol[num], df)
        atomtypes_stringname = Symbol(string(mol_df.molname[num],"atomtypes"))
        @eval(($curr_mol_name) = ($curr_mol))
        @eval(($atomtypes_stringname) = ($atomtypes_list))
    end
end

function string_as_variable_name(str::AbstractString, var::Any)
    str = Symbol(str)
    return @eval (($str) = ($var))
end