export
    load_sdfile,
    write_sdfile

"""
    load_sdfile(fname::AbstractString, ::Type{T} = Float32) -> System{T}

Read an SDFile.
"""
function load_sdfile(fname::AbstractString, ::Type{T} = Float32) where {T <: Real}
    sys = System{T}(basename(fname))
    for mg_mol in MolecularGraph.sdfilereader(fname)
        convert(Molecule{T}, mg_mol; system=sys)
    end
    sys
end

"""
    write_sdfile(fname::AbstractString, ac::AbstractAtomContainer)

Save a 2D projection of an atom container as SDFile.
"""
@inline function write_sdfile(fname::AbstractString, mol::AbstractAtomContainer)
    MolecularGraph.sdfilewriter(fname, [convert(MolecularGraph.SDFMolGraph, mol)])
end

@inline function write_sdfile(fname::AbstractString, sys::System)
    MolecularGraph.sdfilewriter(fname, convert.(MolecularGraph.SDFMolGraph, molecules(sys)))
end
