export System, default_system, parent_system

const _SystemAtomTuple{T} = NamedTuple{
    (fieldnames(AtomTuple{T})...,
        :frame_id, :molecule_id, :chain_id, :fragment_id, :nucleotide_id, :residue_id),
    Tuple{fieldtypes(AtomTuple{T})...,
        Int, MaybeInt, MaybeInt, MaybeInt, MaybeInt, MaybeInt}
}

const _SystemChainTuple = NamedTuple{
    (fieldnames(ChainTuple)..., :molecule_id),
    Tuple{fieldtypes(ChainTuple)..., Int}
}

const _SystemFragmentTuple = NamedTuple{
    (fieldnames(FragmentTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(FragmentTuple)..., Int, Int}
}

const _SystemNucleotideTuple = NamedTuple{
    (fieldnames(NucleotideTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(NucleotideTuple)..., Int, Int}
}

const _SystemResidueTuple = NamedTuple{
    (fieldnames(ResidueTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(ResidueTuple)..., Int, Int}
}

mutable struct System{T} <: AbstractAtomContainer{T}
    name::String
    atoms::DataFrame
    bonds::DataFrame
    molecules::DataFrame
    chains::DataFrame
    fragments::DataFrame
    nucleotides::DataFrame
    residues::DataFrame
    properties::Properties

    _curr_idx::Int

    function System{T}(name::String = "") where T
        new(
            name,
            DataFrame(_SystemAtomTuple{T}[]),
            DataFrame(BondTuple[]),
            DataFrame(MoleculeTuple[]),
            DataFrame(_SystemChainTuple[]),
            DataFrame(_SystemFragmentTuple[]),
            DataFrame(_SystemNucleotideTuple[]),
            DataFrame(_SystemResidueTuple[]),
            Properties(), 
            0
        )
    end
end

System(name::String = "") = System{Float32}(name)

const _default_system = System("default")

@inline function default_system()
    _default_system
end

@inline function _next_idx(sys::System{T}) where T
    sys._curr_idx += 1
end

Base.show(io::IO, ::MIME"text/plain", sys::System) = show(io, sys)
Base.show(io::IO, sys::System) = print(io, 
    "System with ", natoms(sys), " atoms", isempty(sys.name) ? "" : " ($(sys.name))")

#=
    System-aware DataFrame wrapper
=#
struct SystemDataFrame{T} <: AbstractDataFrame
    sytem::System{T}
    df::SubDataFrame
end

DataFrames.describe(sys::SystemDataFrame) = describe(getfield(sys, :df))
DataFrames.index(sys::SystemDataFrame) = DataFrames.index(getfield(sys, :df))
DataFrames.nrow(sys::SystemDataFrame) = nrow(getfield(sys, :df))
DataFrames._check_consistency(sys::SystemDataFrame) = DataFrames._check_consistency(getfield(sys, :df))
Base.getindex(sys::SystemDataFrame, row, col) = getindex(getfield(sys, :df), row, col)
Base.getindex(sys::SystemDataFrame, row::Integer, col::Colon) = getindex(getfield(sys, :df), row, col)
Base.setindex!(sys::SystemDataFrame, val, idx) = setindex!(getfield(sys, :df), val, idx)
Base.setindex!(sys::SystemDataFrame, val, row, col) = setindex!(getfield(sys, :df), val, row, col)
Base.getproperty(sys::SystemDataFrame, name::Symbol) = getproperty(getfield(sys, :df), name)
Base.getproperty(sys::SystemDataFrame, name::Symbol, order::Symbol) = getproperty(getfield(sys, :df), name, order)
Base.setproperty!(sys::SystemDataFrame, col, v) = setproperty!(getfield(sys, :df), col, v)
Base.show(io::IO, ::MIME"text/plain", sys::SystemDataFrame) = show(io, getfield(sys, :df))
Base.show(io::IO, sys::SystemDataFrame) = show(io, getfield(sys, :df))
