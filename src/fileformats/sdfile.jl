export
    load_sdfile,
    write_sdfile

"""
    load_sdfile(fname::String, ::Type{T} = Float32) -> System{T}

Read an SDFile.
"""
function load_sdfile(fname::String, ::Type{T} = Float32) where {T <: Real}
    mg_mols = [m for m in MolecularGraph.sdfilereader(fname)]

    sys = System{T}(basename(fname))

    for mg_mol in mg_mols
        convert(Molecule{T}, mg_mol; system=sys)
    end

    sys
end

"""
    write_sdfile(fname::String, ac::AbstractAtomContainer)

Save a 2D projection of an atom container as SDFile.
"""
@inline function write_sdfile(fname::String, mol::AbstractAtomContainer)
    @warn "write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane..."
    MolecularGraph.sdfilewriter(fname, [convert(MolecularGraph.SDFMolGraph, mol)])
end

@inline function write_sdfile(fname::String, sys::System)
    @warn "write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane..."
    MolecularGraph.sdfilewriter(fname, convert.(MolecularGraph.SDFMolGraph, molecules(sys)))
end
