export load_sdfile

using MolecularGraph: sdfilereader

function load_sdfile(fname::String, T=Float32)
    mg_mols = [m for m in sdfilereader(fname)]

    [convert(Molecule{T}, m) for m in mg_mols]
end