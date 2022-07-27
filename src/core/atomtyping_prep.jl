
using DataFrames, DelimitedFiles, CSV


function df_maker(mapfile::AbstractString)
    source_file = open(mapfile)
    temp_str = ""
    df = DataFrame([[] for _ = (1:9)], 
    ["ATD", "type_name", "residue_names", "atomic_number", "num._attached_atoms", "num._attached_H", "electron_withdrawal_groups", "atomic_property", "CES"])
    
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
