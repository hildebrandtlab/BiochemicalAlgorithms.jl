export Bond, bond_by_idx, bonds, bonds_df, eachbond, nbonds, get_partner, is_bound_to

"""
    $(TYPEDEF)

Mutable representation of an individual bond in a system.

# Fields
 - `idx::Int`
 - `a1::Int`
 - `a2::Int`
 - `order::BondOrderType`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Bond(
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Bond{Float32}` in the default system.

```julia
Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Bond{T}` in the given system.
"""
@auto_hash_equals struct Bond{T} <: AbstractSystemComponent{T}
    _sys::System{T}
    _row::DataFrameRow
end

function Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType, 
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    idx = _next_idx(sys)
    push!(sys._bonds, (idx, a1, a2, order, properties, flags))
    bond_by_idx(sys, idx)
end

@inline function Bond(
    a1::Int,
    a2::Int,
    order::BondOrderType,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    Bond(default_system(), a1, a2, order, properties, flags)
end

function Base.getproperty(bond::Bond, name::Symbol)
    in(name, fieldnames(BondTuple)) && return getproperty(getfield(bond, :_row), name)
    getfield(bond, name)
end

function Base.setproperty!(bond::Bond, name::Symbol, val)
    in(name, fieldnames(BondTuple)) && return setproperty!(getfield(bond, :_row), name, val)
    setproperty!(bond, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", bond::Bond) = show(io, getfield(bond, :_row))
@inline Base.show(io::IO, bond::Bond) = show(io, getfield(bond, :_row))

@inline Base.parent(bond::Bond) = bond._sys
@inline parent_system(bond::Bond) = parent(bond)
# TODO other parent_ functions

"""
    $(TYPEDSIGNATURES)

Returns the `Bond{T}` associated with the given `idx` in `sys`. Returns `nothing` if no such bond
exists.
"""
@inline function bond_by_idx(sys::System{T}, idx::Int) where T
    rn = _row_by_idx(sys._bonds, idx)
    isnothing(rn) ? nothing : Bond{T}(sys, DataFrameRow(sys._bonds, rn, :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's bonds associated with at least one atom
matching the given criteria (value or `missing`). Fields given as `nothing` are ignored. See
[`_atoms`](@ref).
"""
function _bonds(sys::System; kwargs...)
    # FIXME this implementation currently ignores bonds with _two_ invalid atom IDs
    aidx = _atoms(sys; kwargs...).idx
    @rsubset(
        sys._bonds, :a1 in aidx || :a2 in aidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
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

Returns a `SubDataFrame` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
 - `frame_id::Union{Nothing, Int} = 1`: \
Any value other than `nothing` also limits the result to bonds where at least on atom matches \
this frame ID.
"""
@inline function bonds_df(sys::System; kwargs...)
    _bonds(sys; kwargs...)
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
    push!(ac::AbstractAtomContainer, bond::BondTuple)

Creates a new bond in the system associated with the given atom container, based on the given tuple.
The new bond is automatically assigned a new `idx`.
"""
function Base.push!(sys::System{T}, bond::BondTuple) where T
    push!(sys._bonds, (; bond..., idx = _next_idx(sys)))
    sys
end

function get_partner(bond, atom)
    if bond.a1 == atom.idx
        return atom_by_idx(atom._sys, bond.a2)
    elseif bond.a2 == atom.idx
        return atom_by_idx(atom._sys, bond.a1)
    else
        return nothing
    end
end

function is_bound_to(a1::Atom, a2::Atom)
    s = a1._sys

    if s != a2._sys
        return false
    end

    return !isnothing(
        findfirst(
            b -> 
                ((b.a1 == a1.idx) && (b.a2 == a2.idx)) ||
                ((b.a1 == a2.idx) && (b.a2 == a1.idx)), 
            bonds(s)
        )
    )
end
