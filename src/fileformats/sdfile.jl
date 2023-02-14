export load_sdfile, write_sdfile

using MolecularGraph: sdfilereader, sdfilewriter,
    GraphMol, SDFileAtom, SDFileBond

function load_sdfile(fname::String, T=Float32)
    mg_mols = [m for m in sdfilereader(fname)]

    [convert(Molecule{T}, m) for m in mg_mols]
end

function write_sdfile(fname::String, mol::AbstractMolecule)
    mg_mol = convert(GraphMol{SDFileAtom, SDFileBond}, mol)

    sdfilewriter(fname, [mg_mol])
end

function write_sdfile(fname::String, mols::AbstractArray{M}) where {M<:AbstractMolecule}
    mg_mols = [convert(GraphMol{SDFileAtom, SDFileBond}, m) for m in mols]

    sdfilewriter(fname, mg_mols)
end