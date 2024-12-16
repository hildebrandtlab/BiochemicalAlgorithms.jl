export
    BALLIniFile,
    read_ball_ini_file

using CSV
using DataStructures: OrderedDict

@auto_hash_equals mutable struct BALLIniFileSection
    name::String
    properties::OrderedDict{String, String}
    data::DataFrame

    function BALLIniFileSection(name="", properties=OrderedDict{String, String}(), data=DataFrame())
        new(name, properties, data)
    end
end

@auto_hash_equals struct BALLIniFile
    sections::OrderedDict{String, BALLIniFileSection}

    function BALLIniFile(sections=OrderedDict{String, BALLIniFileSection}())
        new(sections)
    end
end

"""
    read_ball_ini_file(fname::String, ::Type{T} = Float32) -> BALLIniFile

Read a file in BALL's old Ini file format and return it as a BALLIniFile.

# Supported keyword arguments
 - `cleanup_keys::Bool = true` simplifies colon-separated key names
   (e.g., `ver:version` becomes `version`, `key:I` becomes `I`, etc.)
"""
function read_ball_ini_file(fname::String, ::Type{T} = Float32; cleanup_keys::Bool = true) where {T <: Real}
    result = BALLIniFile()

    open(fname) do paramfile
        # filter out comments
        params_string = join(filter(l -> !startswith(strip(l), ";"), readlines(paramfile)), "\n")

        # find each section
        for s_match in eachmatch(r"^\[(.*)\]"m, params_string)
            section = BALLIniFileSection()

            s_name = s_match.captures[1]

            section.name = s_name

            # start at the next line
            from = match(r"\n", params_string, s_match.offset)

            if isnothing(from)
                result.sections[s_name] = section
                continue
            end

            # and read until the next empty line
            to = match(r"^\s*$"m, params_string, from.offset+1)

            data = params_string[from.offset+1:(isnothing(to) ? end : to.offset-1)]

            # does the first line contain a header?
            has_header = !isnothing(match(r"^.*:.*$"m, data))

            # are there any property lines?
            all_property_matches = eachmatch(r"^\s*@(.*)=(.*)$"m, data)

            section.properties = OrderedDict{String, String}(
                m[1] => m[2] for m in all_property_matches
            )

            # filter out property lines
            data = join(filter(l -> !startswith(strip(l), "@"), split(data, "\n")), "\n")

            s_df = CSV.read(
                IOBuffer(data), DataFrame;
                header=has_header,
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

            if cleanup_keys
                rename!(
                    s_df,
                    Dict( n => occursin(":", n) ? split(n, ":")[2] : n for n in names(s_df) )
                )
            end

            section.data = s_df

            result.sections[s_name] = section
        end
    end

    result
end
