export AbstractForceFieldParameters, extract_section

abstract type AbstractForceFieldParameters end

function extract_section(params::AbstractForceFieldParameters, section::String)
    if !haskey(params.sections, section)
        @error "extract_section(): Could not extract section $(section) from $(params.filename)!"
    end

    params.sections[section]
end