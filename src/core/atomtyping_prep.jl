using CSV
using DataFrames
using DelimitedFiles

export DEFtoCSV_converter

function DEFtoCSV_converter(mapfile::AbstractString)
    source_file = open(mapfile)
    df = DataFrame([[] for _ = (1:9)],
        ["ATD" , "type_name", "residue_names", "atomic_number", "num_neighbors",
        "num_H_bonds", "electron_withdrawal_groups", "atomic_property", "CES"])
    # [AbstractString, AbstractString, Char, Int64, Int64, Int64, Int64, AbstractString, AbstractString],
    ### possibly add type casting for columns ??
    
    for line in readlines(source_file)
        if lastindex(line)>=3 && line[1:3] == "ATD"
            col_data = lines_for_df(line)
            push!(df, col_data)
        end
    end
    
    out_str = lowercase(mapfile[26:lastindex(mapfile)-4])
    CSV.write("$out_str.csv", df)
    println(df)
    # return df
end

function lines_for_df(line::AbstractString)
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
    for i = (1:9)
        if i >= 4 && i <= 7 && (col_data[i] == "" || col_data[i] == "*")
            col_data[i] = "-1"
        elseif col_data[i] == ""
            col_data[i] = "*"
        end   
    end         
    return col_data
end


