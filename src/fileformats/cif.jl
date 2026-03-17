# CIF 1.1 / CIF 2.0 Parser in Julia
# Includes in-memory model for structured access

@enum CIFVersion v1_1 v2_0
@enum ParserState begin
    Start
    InDataBlock
    InLoopHeader
    InLoopBody
    InMultilineText
    InTripleQuotedString
    ExpectValue
    Error
end

struct CIFLoop
    tags::Vector{String}
    rows::Vector{Vector{String}}
end

mutable struct CIFDataBlock
    name::String
    tags::Dict{String, String}
    loops::Vector{CIFLoop}
end

mutable struct CIFFile
    version::CIFVersion
    blocks::Dict{String, CIFDataBlock}
end

mutable struct CIFParser
    version::CIFVersion
    state::ParserState
    current_tag::Union{Nothing, String}
    loop_tags::Vector{String}
    loop_data::Vector{String}  # flat accumulation of all loop values
    line_number::Int
    in_text::Bool
    text_buffer::Vector{String}
    triple_quote_delim::Union{Nothing, String}
    current_block::Union{Nothing, CIFDataBlock}
    file::CIFFile
    pre_text_state::ParserState  # state before entering text block
end

function CIFParser()
    CIFParser(v1_1, Start, nothing, String[], String[], 0, false, String[], nothing, nothing, CIFFile(v1_1, Dict()), Start)
end

function parse_line!(parser::CIFParser, line::String)
    parser.line_number += 1
    stripped = strip(line)

    if parser.line_number == 1 && startswith(stripped, "#\\#CIF_2.0")
        parser.version = v2_0
        parser.file.version = v2_0
        return
    end

    if isempty(stripped) || startswith(stripped, '#')
        return
    end

    if parser.version == v2_0 && !isnothing(parser.triple_quote_delim)
        if occursin(parser.triple_quote_delim, line)
            push!(parser.text_buffer, split(line, parser.triple_quote_delim)[1])
            store_key_value(parser, parser.current_tag, join(parser.text_buffer, "\n"))
            parser.triple_quote_delim = nothing
            parser.current_tag = nothing
            parser.state = InDataBlock
        else
            push!(parser.text_buffer, line)
        end
        return
    end

    if parser.in_text
        if startswith(line, ";")
            text_value = join(parser.text_buffer, "\n")
            if parser.pre_text_state == InLoopBody
                # Multiline text inside a loop — append to loop data
                push!(parser.loop_data, text_value)
                parser.state = InLoopBody
            else
                store_key_value(parser, parser.current_tag, text_value)
                parser.state = InDataBlock
            end
            parser.current_tag = nothing
            parser.text_buffer = String[]
            parser.in_text = false
        else
            push!(parser.text_buffer, line)
        end
        return
    end

    if startswith(stripped, "data_")
        # Finalize any pending loop before switching data blocks
        finalize_loop!(parser)
        name = stripped[6:end]
        parser.current_block = CIFDataBlock(name, Dict(), CIFLoop[])
        parser.file.blocks[name] = parser.current_block
        parser.state = InDataBlock
        return
    elseif startswith(stripped, "loop_")
        # Finalize any pending loop before starting a new one
        finalize_loop!(parser)
        parser.state = InLoopHeader
        empty!(parser.loop_tags)
        empty!(parser.loop_data)
        return
    elseif parser.state == InLoopHeader && startswith(stripped, "_")
        push!(parser.loop_tags, String(stripped))
        return
    elseif parser.state in (InLoopHeader, InLoopBody) && !startswith(stripped, "_")
        # Check for semicolon text block in loop data
        if startswith(line, ";")
            parser.pre_text_state = InLoopBody
            parser.in_text = true
            parser.text_buffer = [line[2:end]]
            parser.state = InMultilineText
            return
        end
        parser.state = InLoopBody
        values = parse_compound_values(String(stripped))
        append!(parser.loop_data, values)
        return
    elseif startswith(stripped, "_")
        # If we were in a loop, finalize it first
        if parser.state == InLoopBody
            finalize_loop!(parser)
        end
        # Tag-value pair: may be on one line or split across two
        tokens = parse_compound_values(String(stripped))
        if length(tokens) >= 2
            # Inline tag-value: _tag value
            store_key_value(parser, tokens[1], join(tokens[2:end], " "))
            parser.state = InDataBlock
        else
            parser.current_tag = tokens[1]
            parser.state = ExpectValue
        end
        return
    elseif parser.state == ExpectValue
        if parser.version == v2_0 && (startswith(line, "\"\"\"") || startswith(line, "'''"))
            delim = startswith(line, "\"\"\"") ? "\"\"\"" : "'''"
            content = replace(line, delim => "")
            parser.triple_quote_delim = delim
            parser.text_buffer = [content]
            parser.state = InTripleQuotedString
            return
        elseif startswith(line, ";")
            parser.pre_text_state = ExpectValue
            parser.in_text = true
            parser.text_buffer = [line[2:end]]
            parser.state = InMultilineText
            return
        else
            store_key_value(parser, parser.current_tag, String(stripped))
            parser.current_tag = nothing
            parser.state = InDataBlock
            return
        end
    else
        parser.state = Error
        @warn "Unexpected line at $(parser.line_number): $line"
    end
end

function parse_compound_values(line::String)
    result = String[]
    buffer = IOBuffer()
    in_list = false
    in_table = false
    in_quote = false
    quote_char = ' '

    for c in line
        if c in ('"', '\'')
            if in_quote && c == quote_char
                in_quote = false
                push!(result, String(take!(buffer)))
            elseif !in_quote
                in_quote = true
                quote_char = c
            else
                write(buffer, c)
            end
        elseif c == '['
            in_list = true
            write(buffer, c)
        elseif c == ']'
            in_list = false
            write(buffer, c)
            push!(result, String(take!(buffer)))
        elseif c == '{'
            in_table = true
            write(buffer, c)
        elseif c == '}'
            in_table = false
            write(buffer, c)
            push!(result, String(take!(buffer)))
        elseif isspace(c) && !in_quote && !in_list && !in_table
            if buffer.size > 0
                push!(result, String(take!(buffer)))
            end
        else
            write(buffer, c)
        end
    end
    if buffer.size > 0
        push!(result, String(take!(buffer)))
    end
    return result
end

function store_key_value(parser::CIFParser, tag::Union{Nothing, String}, value::String)
    tag === nothing && return
    if parser.current_block !== nothing
        parser.current_block.tags[tag] = value
    end
end

function finalize_loop!(parser::CIFParser)
    if parser.current_block !== nothing && !isempty(parser.loop_tags) && !isempty(parser.loop_data)
        ntags = length(parser.loop_tags)
        nvals = length(parser.loop_data)

        if nvals % ntags != 0
            @warn "CIF loop has $(nvals) values for $(ntags) tags — not evenly divisible"
        end

        # Reshape flat values into rows of ntags columns
        rows = Vector{Vector{String}}()
        for i in 1:ntags:nvals
            last_idx = min(i + ntags - 1, nvals)
            push!(rows, parser.loop_data[i:last_idx])
        end

        push!(parser.current_block.loops, CIFLoop(copy(parser.loop_tags), rows))
    end
    empty!(parser.loop_tags)
    empty!(parser.loop_data)
end

function parse_cif(io::IO)
    parser = CIFParser()
    for line in eachline(io)
        parse_line!(parser, line)
    end
    finalize_loop!(parser)
    return parser.file
end

function parse_cif_file(path::String)
    open(path, "r") do io
        parse_cif(io)
    end
end
