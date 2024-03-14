export
    Bond,
    BondTable,
    bond_by_idx,
    bonds,
    nbonds,
    get_partner,
    get_partners,
    bond_length,
    non_hydrogen_bonds,
    hydrogen_bonds

"""
    $(TYPEDEF)

Tables.jl-compatible representation of system bonds (or a subset thereof). Bond tables can be
generated using [`bonds`](@ref) or filtered from other bond tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `a1::AbstractVector{Int}`
 - `a2::AbstractVector{Int}`
 - `order::AbstractVector{BondOrderType}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
"""
@auto_hash_equals struct BondTable{T} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _bonds(bt::BondTable) = getfield(getfield(bt, :_sys), :_bonds)

@inline function Tables.getcolumn(bt::BondTable, nm::Symbol)
    col = Tables.getcolumn(_bonds(bt), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _bonds(bt)._idx_map[idx], bt._idx)
    )
end

@inline Tables.columnnames(bt::BondTable) = Tables.columnnames(_bonds(bt))
@inline Tables.schema(bt::BondTable) = Tables.schema(_bonds(bt))

@inline function Base.getproperty(bt::BondTable, nm::Symbol)
    hasfield(typeof(bt), nm) && return getfield(bt, nm)
    Tables.getcolumn(bt, nm)
end

@inline function Base.setproperty!(bt::BondTable, nm::Symbol, val)
    if nm in _bond_table_cols_priv || nm in _bond_table_cols_set
        error("BondTable columns cannot be set directly! Did you mean to use broadcast assignment (.=)?")
    end
    if !hasfield(typeof(bt), nm)
        error("type BondTable has no field $nm")
    end
    setfield!(bt, nm, val)
end

@inline function _filter_bonds(f::Function, sys::System)
    BondTable(sys, collect(Int, _filter_select(f, sys._bonds, :idx)))
end

@inline function Base.filter(f::Function, bt::BondTable)
    BondTable(bt._sys, collect(Int, _filter_select(f, bt, :idx)))
end

@inline function Base.iterate(bt::BondTable, st = 1)
    st > length(bt) ?
        nothing :
        (bond_by_idx(bt._sys, bt._idx[st]), st + 1)
end

@inline Base.eltype(::BondTable{T}) where T = Bond{T}
@inline Base.size(bt::BondTable) = (length(bt._idx), length(Tables.columnnames(bt)))
@inline Base.getindex(bt::BondTable, i::Int) = bond_by_idx(bt._sys, bt._idx[i])
@inline Base.getindex(bt::BondTable, ::Colon) = bt

@inline function Base.getindex(bt::BondTable, I)
    BondTable(bt._sys, collect(Int, map(i -> bt._idx[i], I)))
end

"""
    $(TYPEDEF)

Mutable representation of an individual bond in a system.

# Public fields
 - `idx::Int`
 - `a1::Int`
 - `a2::Int`
 - `order::BondOrderType`

# Private fields
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType;
    # keyword arguments
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Bond{T}` in the given system.

```julia
Bond(
    a1::Int,
    a2::Int,
    order::BondOrderType;
    # keyword arguments
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Bond{Float32}` in the default system.
"""
@auto_hash_equals struct Bond{T} <: AbstractSystemComponent{T}
    _sys::System{T}
    _row::_BondTableRow
end

@inline function Bond(
    sys::System{T}, 
    a1::Int, 
    a2::Int, 
    order::BondOrderType;
    kwargs...
) where T
    idx = _next_idx(sys)
    push!(sys._bonds, idx, a1, a2, order; kwargs...)
    bond_by_idx(sys, idx)
end

@inline function Bond(
    a1::Int,
    a2::Int,
    order::BondOrderType;
    kwargs...
)
    Bond(default_system(), a1, a2, order; kwargs...)
end

@inline function Bond(
    ac::AbstractAtomContainer,
    a1::Int,
    a2::Int,
    order::BondOrderType;
    kwargs...
)
    Bond(parent(ac), a1, a2, order; kwargs...)
end

@inline function Base.getproperty(bond::Bond, name::Symbol)
    hasfield(typeof(bond), name) && return getfield(bond, name)
    getproperty(getfield(bond, :_row), name)
end

@inline function Base.setproperty!(bond::Bond, name::Symbol, val)
    hasfield(typeof(bond), name) && return setfield!(bond, name, val)
    setproperty!(getfield(bond, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", bond::Bond) = show(io, bond)
@inline function Base.show(io::IO, bond::Bond)
    print(io, "$(typeof(bond)): ")
    show(io, NamedTuple(bond._row))
end

@inline Base.parent(bond::Bond) = bond._sys
@inline parent_system(bond::Bond) = parent(bond)
# TODO other parent_ functions

"""
    $(TYPEDSIGNATURES)

Returns the `Bond{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
bond exists.
"""
@inline function bond_by_idx(sys::System{T}, idx::Int) where T
    Bond{T}(sys, _row_by_idx(sys._bonds, idx))
end

"""
    bonds(::Chain)
    bonds(::Fragment)
    bonds(::Molecule)
    bonds(::Nucleotide)
    bonds(::Residue)
    bonds(::System)

Returns a `BondTable{T}` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
function bonds(sys::System; kwargs...)
    # FIXME this implementation currently ignores bonds with _two_ invalid atom IDs
    aidx = Set(atoms(sys; kwargs...).idx)
    _filter_bonds(
        bond -> bond.a1 in aidx || bond.a2 in aidx,
        sys
    )
end

"""
    nbonds(::Chain)
    nbonds(::Fragment)
    nbonds(::Molecule)
    nbonds(::Nucleotide)
    nbonds(::Residue)
    nbonds(::System)

Returns the number of bonds in the given atom container where at least one associated atom
is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
function nbonds(sys::System{T}; kwargs...) where T
    length(bonds(sys; kwargs...))
end

"""
    push!(::AbstractAtomContainer, ::Bond{T})

Creates a copy of the given bond in the system associated with the given atom container.
The new bond is automatically assigned a new `idx`.
"""
@inline function Base.push!(sys::System{T}, bond::Bond{T}) where T
    Bond(sys, bond.a1, bond.a2, bond.order;
        properties = bond.properties,
        flags = bond.flags
    )
    sys
end

@inline function Base.push!(ac::AbstractAtomContainer{T}, bond::Bond{T}) where T
    push!(parent(ac), bond)
    ac
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

function get_partners(bond)
    s = parent(bond)
    atom_by_idx(s, bond.a1), atom_by_idx(s, bond.a2)
end

@inline function bond_length(bond)
    s = parent(bond)
    distance(atom_by_idx(s, bond.a1), atom_by_idx(s, bond.a2))
end

@inline non_hydrogen_bonds(ac) = filter(b -> !haskey(b.properties, :TYPE__HYDROGEN), bonds(ac))
@inline hydrogen_bonds(ac) = filter(b -> has_key(b.properties, :TYPE__HYDROGEN), bonds(ac))
