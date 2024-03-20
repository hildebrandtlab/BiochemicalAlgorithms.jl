export
    load_sdfile,
    write_sdfile

function load_sdfile(fname::String, T=Float32)
    mg_mols = [m for m in MolecularGraph.sdfilereader(fname)]

    sys = System{T}(fname)

    for mg_mol in mg_mols
        convert(Molecule{T}, mg_mol; system=sys)
    end

    sys
end

@inline function write_sdfile(fname::String, mol::AbstractAtomContainer)
    @warn "write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane..."
    MolecularGraph.sdfilewriter(fname, [convert(MolecularGraph.SDFMolGraph, mol)])
end

@inline function write_sdfile(fname::String, sys::System)
    @warn "write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane..."
    MolecularGraph.sdfilewriter(fname, convert.(MolecularGraph.SDFMolGraph, molecules(sys)))
end
