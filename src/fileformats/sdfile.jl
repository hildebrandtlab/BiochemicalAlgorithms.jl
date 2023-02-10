export load_sdfile

using MolecularGraph: sdfilereader

function load_sdfile(fname::String, T=Float32)
    mg_mols = [m for m in sdfilereader(fname)]

    sys = System{T}(fname)

    for mg_mol in mg_mols
        convert(Molecule{T}, mg_mol; system=sys)
    end

    sys
end