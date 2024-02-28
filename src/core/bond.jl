export
    Bond,
    BondTable,
    bond_by_idx,
    bonds,
    bonds_df,
    nbonds,
    get_partner,
    get_partners,
    bond_length,
    non_hydrogen_bonds,
    hydrogen_bonds,
    non_hydrogen_bonds_df,
    hydrogen_bonds_df

@auto_hash_equals struct BondTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _bonds(bt::BondTable) = getproperty(getfield(bt, :_sys), :_bonds)

@inline Tables.istable(::Type{<: BondTable}) = true
@inline Tables.columnaccess(::Type{<: BondTable}) = true
@inline Tables.columns(bt::BondTable) = bt

@inline function Tables.getcolumn(bt::BondTable, nm::Symbol)
    col = Tables.getcolumn(_bonds(bt), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _bonds(bt)._idx_map[idx], getfield(bt, :_idx))
    )
end

@inline function Base.getproperty(bt::BondTable, nm::Symbol)
    hasfield(typeof(bt), nm) && return getfield(bt, nm)
    Tables.getcolumn(bt, nm)
end

@inline Tables.getcolumn(bt::BondTable, i::Int) = Tables.getcolumn(bt, Tables.columnnames(bt)[i])
@inline Tables.columnnames(bt::BondTable) = Tables.columnnames(_bonds(bt))
@inline Tables.schema(bt::BondTable) = Tables.schema(_bonds(bt))

@inline Base.size(bt::BondTable) = (length(getfield(bt, :_idx)), length(_bond_table_cols))
@inline Base.size(bt::BondTable, dim) = size(bt)[dim]
@inline Base.length(bt::BondTable) = size(bt, 1)

function Base.push!(bt::BondTable, t::BondTuple; kwargs...)
    sys = getfield(bt, :_sys)
    push!(sys, t; kwargs...)
    push!(getfield(bt, :_idx), sys._curr_idx)
    bt
end

@inline function _filter_bonds(f::Function, sys::System{T}) where T
    BondTable(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._bonds),
        :idx
    )))
end

@inline function Base.filter(f::Function, bt::BondTable)
    BondTable(getfield(bt, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, bt),
        :idx
    )))
end

@inline function Base.iterate(bt::BondTable, st = 1)
    st > length(bt) ?
        nothing :
        (bond_by_idx(getfield(bt, :_sys), getfield(bt, :_idx)[st]), st + 1)
end
@inline Base.eltype(::BondTable{T}) where T = Bond{T}
@inline Base.getindex(bt::BondTable{T}, i::Int) where T = bond_by_idx(getfield(bt, :_sys), getfield(bt, :_idx)[i])
@inline Base.keys(bt::BondTable) = LinearIndices((length(bt),))

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
    _row::_BondTableRow
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
    push!(sys._bonds, BondTuple(a1, a2, order;
        idx = idx,
        properties = properties,
        flags = flags
    ))
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

@inline Tables.rows(bt::BondTable) = bt
@inline Tables.getcolumn(bond::Bond, name::Symbol) = Tables.getcolumn(getfield(bond, :_row), name)

@inline function Base.getproperty(bond::Bond, name::Symbol)
    hasfield(typeof(bond), name) && return getfield(bond, name)
    getproperty(getfield(bond, :_row), name)
end

@inline function Base.setproperty!(bond::Bond, name::Symbol, val)
    hasfield(typeof(bond), name) && return setfield!(bond, name, val)
    setproperty!(getfield(bond, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", bond::Bond) = show(io, getfield(bond, :_row))
@inline Base.show(io::IO, bond::Bond) = show(io, getfield(bond, :_row))

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
    bonds(::Protein)
    bonds(::Residue)
    bonds(::System)

Returns a `BondTable{T}` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
function bonds(sys::System; kwargs...)
    # FIXME this implementation currently ignores bonds with _two_ invalid atom IDs
    aidx = Set(_filter_select(atoms(sys; kwargs...), :idx))
    _filter_bonds(
        bond -> bond.a1 in aidx || bond.a2 in aidx,
        sys
    )
end

"""
    bonds_df(::Chain)
    bonds_df(::Fragment)
    bonds_df(::Molecule)
    bonds_df(::Nucleotide)
    bonds_df(::Protein)
    bonds_df(::Residue)
    bonds_df(::System)

Returns a `DataFrame` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
@inline function bonds_df(sys::System; kwargs...)
    DataFrame(bonds(sys; kwargs...))
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
See [`atoms`](@ref)
"""
function nbonds(sys::System{T}; kwargs...) where T
    length(bonds(sys; kwargs...))
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

@inline non_hydrogen_bonds_df(ac) = filter(b->:TYPE__HYDROGEN ∉ keys(b.properties), bonds_df(ac))
@inline hydrogen_bonds_df(ac) = filter(b->:TYPE__HYDROGEN ∈ keys(b.properties), bonds_df(ac))
