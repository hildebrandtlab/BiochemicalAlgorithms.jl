using CSV
using DataFrames
using DelimitedFiles

export df_maker


function df_maker(mapfile::AbstractString)
    source_file = open(mapfile)
    df = DataFrame([[] for _ = (1:9)], 
    ["ATD" , "type_name", "residue_names", "atomic_number", "num_neighbors",
     "num_H_bonds", "electron_withdrawal_groups", "atomic_property", "CES"])
    ### possibly add type casting for columns ??
    
    for line in readlines(source_file)
        if lastindex(line)>=3 && line[1:3] == "ATD"
            col_data = column_info_from_lines(line)
            push!(df, col_data)
        end
    end
    CSV.write("atomtype_gff.csv", df)
    return df
end

function column_info_from_lines(line::AbstractString)
    col_data = ["" for _ = (1:9)]
    current_col = 1
    for i = (1:lastindex(line))
        if line[i] == &
            break
        elseif line[i] != ' ' && line[i] != '\t' && line[i] != '&'
            col_data[current_col] = string(col_data[current_col], line[i])
        elseif i >= 3 && (line[i] == ' ' || line[i] == '\t') && line[i-1] != ' ' && line[i-1] != '\t'
            current_col += 1
        end
    end
    return col_data
end


