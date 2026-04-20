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
    Bond{T} <: AbstractAtomContainer{T}

Mutable representation of an individual bond in a system.

# Public fields
 - `idx::Int`
 - `order::BondOrderType`

# Private fields
 - `properties::Properties`
 - `flags::Flags`
 - `atom1_idx::Int`
 - `atom2_idx::Int`

# Constructors
```julia
Bond(
    a1::Atom{T},
    a2::Atom{T},
    order::BondOrderType;
    # keyword arguments
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Bond{T}` for the given atoms. Both atoms must belong to the same system.
"""
const Bond{T} = AtomContainer{T, :Bond}

@inline function Bond(
    a1::Atom{T},
    a2::Atom{T},
    order::BondOrderType;
    kwargs...
) where T
    sys = parent(a1)
    @assert sys === parent(a2) "given atoms must belong to the same system to form a bond"
    idx = _next_idx!(sys)
    push!(sys._bonds, idx, order, a1.idx, a2.idx; kwargs...)
    bond_by_idx(sys, idx)
end

"""
    BondTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system bonds (or a subset thereof). Bond tables can be
generated using [`bonds`](@ref) or filtered from other bond tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `order::AbstractVector{BondOrderType}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `atom1_idx::AbstractVector{Int}`
 - `atom2_idx::AbstractVector{Int}`
"""
const BondTable{T} = SystemComponentTable{T, Bond{T}}

@inline function _wrap_bonds(sys::System{T}) where T
    BondTable{T}(sys, getfield(getfield(sys, :_bonds), :idx))
end

@inline function _filter_bonds(f::Function, sys::System{T}) where T
    BondTable{T}(sys, _filter_idx(f, sys._bonds))
end

@inline _table(sys::System{T}, ::Type{Bond{T}}) where T = sys._bonds

@inline function _hascolumn(::Type{<: Bond}, nm::Symbol)
    _hascolumn(_BondTable, nm)
end

"""
    bond_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Bond{T}

Returns the `Bond{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
bond exists.
"""
@inline function bond_by_idx(sys::System{T}, idx::Int) where T
    _rowno_by_idx(_table(sys, Bond{T}), idx) # check idx
    Bond{T}(sys, idx)
end

@inline function bond_by_idx(idx::Int)
    bond_by_idx(default_system(), idx)
end

"""
    bonds(::Chain)
    bonds(::Fragment)
    bonds(::Molecule)
    bonds(::System = default_system())

Returns a `BondTable{T}` containing all bonds of the given atom container where at least one
associated atom is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
function bonds(sys::System = default_system(); kwargs...)
    # FIXME this implementation currently ignores bonds with _two_ invalid atom IDs
    aidx = Set(atoms(sys; kwargs...).idx)
    _filter_bonds(
        bond -> bond.atom1_idx in aidx || bond.atom2_idx in aidx,
        sys
    )
end

@doc raw"""
    bonds(::ChainTable)
    bonds(::FragmentTable)
    bonds(::MoleculeTable)

Returns a `BondTable{T}` containing all bonds where at least one associated atom belongs to
the given table.
""" bonds(::SystemComponentTable)

"""
    nbonds(::Chain)
    nbonds(::Fragment)
    nbonds(::Molecule)
    nbonds(::System = default_system())

Returns the number of bonds in the given atom container where at least one associated atom
is contained in the same container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
@inline function nbonds(sys::System = default_system(); kwargs...)
    length(bonds(sys; kwargs...))
end

"""
    nbonds(::BondTable)
    nbonds(::ChainTable)
    nbonds(::FragmentTable)
    nbonds(::MoleculeTable)

Returns the number of bonds where at least one associated atom belongs to the given table.
"""
@inline function nbonds(bt::BondTable)
    length(bt)
end

"""
    delete!(::Bond)
    delete!(::BondTable)
    delete!(::BondTable, idx::Int)

Removes the given bond(s) from the associated system.
"""
@inline function Base.delete!(bond::Bond)
    delete!(parent(bond)._bonds, bond.idx)
    nothing
end

function Base.delete!(bt::BondTable)
    delete!(_table(bt), bt._idx)
    empty!(bt._idx)
    bt
end

function Base.delete!(bt::BondTable, idx::Int)
    idx in bt._idx || throw(KeyError(idx))
    delete!(bond_by_idx(bt._sys, idx))
    deleteat!(bt._idx, findall(i -> i == idx, bt._idx))
    bt
end

function get_partner(bond, atom)
    if bond.atom1_idx == atom.idx
        return atom_by_idx(atom._sys, bond.atom2_idx)
    elseif bond.atom2_idx == atom.idx
        return atom_by_idx(atom._sys, bond.atom1_idx)
    else
        return nothing
    end
end

function get_partners(bond)
    s = parent(bond)
    atom_by_idx(s, bond.atom1_idx), atom_by_idx(s, bond.atom2_idx)
end

@inline function bond_length(bond)
    s = parent(bond)
    distance(atom_by_idx(s, bond.atom1_idx), atom_by_idx(s, bond.atom2_idx))
end

@inline non_hydrogen_bonds(ac) = filter(b -> !has_flag(b, :TYPE__HYDROGEN), bonds(ac))
@inline hydrogen_bonds(ac) = filter(b -> has_flag(b, :TYPE__HYDROGEN), bonds(ac))
