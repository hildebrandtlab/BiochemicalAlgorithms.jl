using AutoHashEquals
using CSV
export MMFF94Parameters

@auto_hash_equals struct MMFF94Parameters{T<:Real}
    sections::Dict{String, DataFrame}
    radii::Vector{T}
    electronegativities::Vector{T}

    function MMFF94Parameters{T}(
            path::String = ball_data_path("forcefields/MMFF94/mmff94.ini")) where {T<:Real}
        sections = Dict{String, DataFrame}()
        open(path) do paramfile
            params_string = read(paramfile, String)

            # find each section
            for s_match in eachmatch(r"^\[(.*)\]"m, params_string)
                s_name = s_match.captures[1]

                # start at the next line
                from = match(r"\n", params_string, s_match.offset)

                # and read until the next empty line
                to = match(r"^\s*$"m, params_string, from.offset+1)

                data = params_string[from.offset+1:to.offset-1]

                s_df = CSV.read(
                    IOBuffer(data), DataFrame; 
                    header=true, 
                    comment="#", 
                    ignoreemptyrows=true, 
                    delim=" ", 
                    ignorerepeated=true, 
                    silencewarnings=true
                )

                # now, change all Float64? columns to T? columns;
                # CSV's typemap argument should do that for us, but
                # right now, this has a bad performance regression
                if T != Float64
                    for (i,t) in enumerate(eltype.(eachcol(s_df)))
                        if (t == Float64) || (t == Union{Missing, Float64})
                            s_df[!, i] = map(d -> ismissing(d) ? d : T(d), s_df[:, i])
                        end
                    end
                end

                sections[s_name] = s_df
            end
        end

        # see http://www.ccl.net/cca/data/MMFF94/
        radii = T.([
             0.33, 0.0,
             1.34, 0.90, 0.81, 0.77, 0.73, 0.72, 0.74, 0.0,
             1.54, 1.30, 1.22, 1.15, 1.09, 1.03, 1.01, 0.0,
             1.96, 1.74,
             1.44, 1.36, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.38, 1.31,
             1.19, 1.20, 1.20, 1.16, 1.15, 0.0,
             2.11, 1.92,
             1.62, 1.48, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.53, 1.48,
             1.46, 1.40, 1.41, 1.35, 1.33, 0.0
        ])
        
        electronegativities = T.([
             2.20, 0.0,
             0.97, 1.47, 2.01, 2.5, 3.07, 3.5, 4.10, 0.0,
             1.01, 1.23, 1.47, 1.74, 2.06, 2.44, 2.83, 0.0,
             0.91, 1.04,
             1.3, 1.5, 1.6, 1.6, 1.5,
             1.8, 1.8, 1.8, 1.9, 1.6,
             1.82, 2.02, 2.20, 2.48, 2.74, 0.0,
             0.89, 0.99,
             1.3, 1.4, 1.6, 1.8, 1.9,
             2.2, 2.2, 2.2, 1.9, 1.7,
             1.49, 1.72, 1.82, 2.01, 2.21, 0.0
        ])

        new(sections, radii, electronegativities)
    end
end

function MMFF94Parameters(path::String = ball_data_path("forcefields/MMFF94/mmff94.ini"))
    MMFF94Parameters{Float32}(path)
end

Base.show(io::IO, mmff_param::MMFF94Parameters{T}) where {T<:Real} = 
    print(io, 
        "MMFF94 parameters with $(length(mmff_param.sections)) sections.")