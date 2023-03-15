export Bond, bonds, bonds_df, eachbond, nbonds

struct Bond{T}
    sys::System{T}
    row::DataFrameRow
end

function Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties()
) where T
    idx = _next_idx(sys)
    push!(sys.bonds, (idx, a1, a2, order, properties))
    _bond_by_idx(sys, idx)
end

@inline function Bond(a1::Int, a2::Int, order::BondOrderType, properties::Properties = Properties())
    Bond(default_system(), a1, a2, order, properties)
end

function Base.getproperty(bond::Bond, name::Symbol)
    in(name, fieldnames(BondTuple)) && return getproperty(getfield(bond, :row), name)
    getfield(bond, name)
end

function Base.setproperty!(bond::Bond, name::Symbol, val)
    in(name, fieldnames(BondTuple)) && return setproperty!(getfield(bond, :row), name, val)
    setproperty!(bond, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", bond::Bond) = show(io, getfield(bond, :row))
@inline Base.show(io::IO, bond::Bond) = show(io, getfield(bond, :row))

@inline Base.parent(bond::Bond) = bond.sys
@inline parent_system(bond::Bond) = parent(bond)
# TODO other parent_ functions

@inline function _bond_by_idx(sys::System{T}, idx::Int) where T
    Bond{T}(sys, DataFrameRow(sys.bonds, findfirst(sys.bonds.idx .== idx), :))
end

function _bonds(sys::System; kwargs...)
    aidx = _atoms(sys; kwargs...).idx
    filter(row -> (row.a1 in aidx || row.a2 in aidx), sys.bonds, view = true)
end

@inline function bonds(sys::System; kwargs...)
    collect(eachbond(sys; kwargs...))
end

@inline function bonds_df(sys::System; kwargs...)
    SystemDataFrame(sys, _bonds(sys; kwargs...))
end

function eachbond(sys::System{T}; kwargs...) where T
    (Bond{T}(sys, row) for row in eachrow(_bonds(sys; kwargs...)))
end

function nbonds(sys::System{T}; kwargs...) where T
    nrow(_bonds(sys; kwargs...))
end

function Base.push!(sys::System{T}, bond::BondTuple) where T
    push!(sys.bonds, _with_idx(bond, _next_idx(sys)))
    sys
end
