
using DataFrames

export load_atomtyping_DEF

function load_atomtyping_DEF(mapfile::AbstractString)
    df = DataFrame([Vector{String}(), Vector{String}(), 
                    Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), 
                    Vector{String}(), Vector{String}()],
                    ["type_name", "residue_names", "atomic_number", "num_neighbors",
                    "num_H_bonds", "electron_withdrawal_groups", "atomic_property", "CES"])
    
    for line in readlines(mapfile)
        if lastindex(line)>=3 && line[1:3] == "ATD"
            push!(df, lines_for_df(line))
        end
    end
    
    return df
end

function lines_for_df(line::AbstractString)
    col_data_list = Vector{Union{AbstractString, Int}}()
    line_data = split(line)
    for i = (2:9)
        if i <= 3
            if (!isassigned(line_data, i) || line_data[i] == "&")
                append!(col_data_list, ["*"])
            else
                append!(col_data_list, [line_data[i]]) 
            end
        elseif i > 3 && i < 8 
            if isassigned(line_data, i)
                if (line_data[i] == "*" || line_data[i] == "&")
                    append!(col_data_list, [-1])
                else
                    append!(col_data_list, [parse(Int, line_data[i])])
                end
            else
                append!(col_data_list, [-1])
            end
        elseif i >= 8
            if isassigned(line_data, i)
                if (line_data[i] == "*" || line_data[i] == "&")
                    append!(col_data_list, ["*"])
                else
                    append!(col_data_list, [line_data[i]]) 
                end
            else 
                append!(col_data_list, ["*"])
            end            
        end 
    end   
    return col_data_list
end


