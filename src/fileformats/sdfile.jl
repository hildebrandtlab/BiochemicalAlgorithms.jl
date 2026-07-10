export
    load_sdfile,
    write_sdfile

"""
    load_sdfile(io::IO, ::Type{T} = Float32) -> System{T}
    load_sdfile(filename::AbstractString, ::Type{T} = Float32) -> System{T}

Read an SDFile.
"""
function load_sdfile(fname_io::Union{AbstractString, IO}, ::Type{T} = Float32) where {T <: Real}
    sys = System{T}()
    for mg_mol in MolecularGraph.sdfilereader(fname_io)
        convert(Molecule{T}, mg_mol; system=sys)
    end
    sys
end

"""
    write_sdfile(io::IO, ac::AbstractAtomContainer)
    write_sdfile(filename::AbstractString, ac::AbstractAtomContainer)

Save a 2D projection of an atom container as SDFile.
"""
@inline function write_sdfile(fname_io::Union{AbstractString, IO}, mol::AbstractAtomContainer)
    MolecularGraph.sdfilewriter(fname_io, [convert(MolecularGraph.SDFMolGraph, mol)])
end

@inline function write_sdfile(fname_io::Union{AbstractString, IO}, sys::System)
    MolecularGraph.sdfilewriter(fname_io, convert.(MolecularGraph.SDFMolGraph, molecules(sys)))
end
