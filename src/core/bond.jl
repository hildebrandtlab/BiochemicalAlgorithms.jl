export Bond, bonds, bonds_df, eachbond, nbonds

"""
    $(TYPEDEF)

Mutable representation of an individual bond in a system.

# Fields
 - `idx::Int`
 - `a1::Int`
 - `a2::Int`
 - `order::BondOrderType`
 - `properties::Properties`

# Constructors
```julia
Bond(
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties()
)
```
Creates a new `Bond{Float32}` in the default system.

```julia
Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties()
)
```
Creates a new `Bond{T}` in the given system.
"""
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

"""
    $(TYPEDSIGNATURES)

Returns the `Bond{T}` associated with the given `idx` in `sys`.
"""
@inline function _bond_by_idx(sys::System{T}, idx::Int) where T
    Bond{T}(sys, DataFrameRow(sys.bonds, findfirst(sys.bonds.idx .== idx), :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's bonds associated with at least one atom
matching the given criteria (value or `missing`). Fields given as `nothing` are ignored. See
[`_atoms`](@ref).
"""
function _bonds(sys::System; kwargs...)
    aidx = _atoms(sys; kwargs...).idx
    filter(row -> (row.a1 in aidx || row.a2 in aidx), sys.bonds, view = true)
end

"""
bonds(::Chain)
bonds(::Fragment)
bonds(::Molecule)
bonds(::Nucleotide)
bonds(::Protein)
bonds(::Residue)
bonds(::System)

Returns a `Vector{Bond{T}}` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
 - `frame_id::Union{Nothing, Int} = 1`: \
Any value other than `nothing` also limits the result to bonds where at least on atom matches \
this frame ID.
"""
@inline function bonds(sys::System; kwargs...)
    collect(eachbond(sys; kwargs...))
end

"""
bonds_df(::Chain)
bonds_df(::Fragment)
bonds_df(::Molecule)
bonds_df(::Nucleotide)
bonds_df(::Protein)
bonds_df(::Residue)
bonds_df(::System)

Returns a `SystemDataFrame{T}` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
 - `frame_id::Union{Nothing, Int} = 1`: \
Any value other than `nothing` also limits the result to bonds where at least on atom matches \
this frame ID.
"""
@inline function bonds_df(sys::System; kwargs...)
    SystemDataFrame(sys, _bonds(sys; kwargs...))
end

"""
eachbond(::Chain)
eachbond(::Fragment)
eachbond(::Molecule)
eachbond(::Nucleotide)
eachbond(::Protein)
eachbond(::Residue)
eachbond(::System)

Returns a `Bond{T}` generator for all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
 - `frame_id::Union{Nothing, Int} = 1`: \
Any value other than `nothing` also limits the result to bonds where at least on atom matches \
this frame ID.
"""
function eachbond(sys::System{T}; kwargs...) where T
    (Bond{T}(sys, row) for row in eachrow(_bonds(sys; kwargs...)))
end

"""
nbonds(::Chain)
nbonds(::Fragment)
nbonds(::Molecule)
nbonds(::Nucleotide)
nbonds(::Protein)
nbonds(::Residue)
nbonds(::System)

Returns the number of bonds in the given atom container where at least one associated atom
is contained in the same container.

# Supported keyword arguments
 - `frame_id::Union{Nothing, Int} = 1`: \
Any value other than `nothing` also limits the result to bonds where at least on atom matches \
this frame ID.
"""
function nbonds(sys::System{T}; kwargs...) where T
    nrow(_bonds(sys; kwargs...))
end

"""
    push!(::Fragment{T}, atom::BondTuple{T})
    push!(::Molecule{T}, atom::BondTuple{T})
    push!(::Nucleotide{T}, atom::BondTuple{T})
    push!(::Protein{T}, atom::BondTuple{T})
    push!(::Residue{T}, atom::BondTuple{T})
    push!(::System{T}, atom::BondTuple{T})

Creates a new bond in the system associated with the given atom container, based on the given tuple.
The new bond is automatically assigned a new `idx`.
"""
function Base.push!(sys::System{T}, bond::BondTuple) where T
    push!(sys.bonds, _with_idx(bond, _next_idx(sys)))
    sys
end
