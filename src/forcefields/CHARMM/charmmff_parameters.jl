export
    CharmmFFParameters

@auto_hash_equals struct CharmmFFParameters <: AbstractForceFieldParameters
    filename::String
    sections::OrderedDict{String, BALLIniFileSection}

    function CharmmFFParameters(
        filename::AbstractString = ball_data_path("forcefields/CHARMM/param19_eef1.ini"),
        ::Type{T} = Float32
    ) where {T <: Real}
        # the CHARMM parameter files contain empty lines within sections
        parameters_ini = read_ball_ini_file(filename, T; blank_lines_end_sections = false)
        new(filename, parameters_ini.sections)
    end
end
