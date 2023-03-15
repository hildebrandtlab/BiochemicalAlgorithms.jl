export Chain, chains, chains_df, eachchain, nchains, parent_chain

struct Chain{T}
    sys::System{T}
    row::DataFrameRow
end

function Chain(
    mol::Molecule{T},
    name::String = "",
    properties::Properties = Properties()
) where T
    sys = mol.sys
    idx = _next_idx(sys)
    push!(sys.chains, (idx, name, properties, mol.idx))
    _chain_by_idx(sys, idx)
end

function Base.getproperty(chain::Chain, name::Symbol)
    in(name, fieldnames(ChainTuple)) && return getproperty(getfield(chain, :row), name)
    getfield(chain, name)
end

function Base.setproperty!(chain::Chain, name::Symbol, val)
    in(name, fieldnames(ChainTuple)) && return setproperty!(getfield(chain, :row), name, val)
    setfield!(chain, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", chain::Chain) = show(io, getfield(chain, :row))
@inline Base.show(io::IO, chain::Chain) = show(io, getfield(chain, :row))

@inline Base.parent(chain::Chain) = chain.sys
@inline parent_system(chain::Chain) = parent(chain)
@inline parent_molecule(chain::Chain) = _molecule_by_idx(chain.sys, chain.row.molecule_id)

@inline function _chain_by_idx(sys::System{T}, idx::Int) where T
    Chain{T}(sys, DataFrameRow(sys.chains, findfirst(sys.chains.idx .== idx), :))
end

function _chains(sys::System; molecule_id::Union{Nothing, Int} = nothing)
    isnothing(molecule_id) && return sys.chains

    get(
        groupby(sys.chains, :molecule_id),
        (molecule_id = molecule_id,),
        DataFrame(_SystemChainTuple[])
    )
end

@inline function chains(sys::System; kwargs...)
    collect(eachchain(sys; kwargs...))
end

@inline function chains_df(sys::System; kwargs...)
    SystemDataFrame(sys, view(_chains(sys; kwargs...), :, 1:length(fieldnames(ChainTuple))))
end

@inline function eachchain(sys::System{T}; kwargs...) where T
    (Chain{T}(sys, row) for row in eachrow(_chains(sys; kwargs...)))
end

@inline function nchains(sys::System; kwargs...)
    nrow(_chains(sys; kwargs...))
end

#=
    Molecule chains
=#
@inline _chains(mol::Molecule) = _chains(mol.sys; molecule_id = mol.idx)
@inline chains(mol::Molecule) = chains(mol.sys, molecule_id = mol.idx)
@inline chains_df(mol::Molecule) = chains_df(mol.sys, molecule_id = mol.idx)
@inline eachchain(mol::Molecule) = eachchain(mol.sys, molecule_id = mol.idx)
@inline nchains(mol::Molecule) = nchains(mol.sys, molecule_id = mol.idx)

@inline function Base.push!(mol::Molecule, chain::ChainTuple)
    push!(sys.chains, (_with_idx(chain, _next_idx(sys))..., mol.idx))
    mol
end

#=
    Chain atoms
=#
@inline _atoms(chain::Chain; kwargs...) = _atoms(chain.sys; chain_id = chain.idx, kwargs...)
@inline atoms(chain::Chain; kwargs...) = atoms(chain.sys; chain_id = chain.idx, kwargs...)
@inline atoms_df(chain::Chain; kwargs...) = atoms_df(chain.sys; chain_id = chain.idx, kwargs...)
@inline eachatom(chain::Chain; kwargs...) = eachatom(chain.sys; chain_id = chain.idx, kwargs...)
@inline natoms(chain::Chain; kwargs...) = natoms(chain.sys; chain_id = chain.idx, kwargs...)

#=
    Chain bonds
=#
@inline _bonds(chain::Chain; kwargs...) = _bonds(chain.sys; chain_id = chain.idx, kwargs...)
@inline bonds(chain::Chain; kwargs...) = bonds(chain.sys; chain_id = chain.idx, kwargs...)
@inline bonds_df(chain::Chain; kwargs...) = bonds_df(chain.sys; chain_id = chain.idx, kwargs...)
@inline eachbond(chain::Chain; kwargs...) = eachbond(chain.sys; chain_id = chain.idx, kwargs...)
@inline nbonds(chain::Chain; kwargs...) = nbonds(chain.sys; chain_id = chain.idx, kwargs...)
