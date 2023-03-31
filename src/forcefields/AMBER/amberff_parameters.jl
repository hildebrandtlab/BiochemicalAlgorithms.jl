export AmberFFParameters

@auto_hash_equals struct AmberFFParameters <: AbstractForceFieldParameters
    filename::String
    sections::OrderedDict{String, BALLIniFileSection}

    function AmberFFParameters(filename::String = ball_data_path("forcefields/AMBER/amber96.ini"), T=Float32)
        parameters_ini = read_ball_ini_file(filename, T)

        new(filename, parameters_ini.sections)
    end
end

