export
    Atom,
    AtomTable,
    atom_by_idx,
    atom_by_name,
    atoms,
    atoms_df,
    eachatom,
    natoms,
    has_property,
    get_property,
    set_property,
    get_full_name,
    is_bound_to,
    is_geminal,
    is_vicinal

@auto_hash_equals struct AtomTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _atoms(at::AtomTable) = getproperty(getfield(at, :_sys), :_atoms)

@inline Tables.istable(::Type{<: AtomTable}) = true
@inline Tables.columnaccess(::Type{<: AtomTable}) = true
@inline Tables.columns(at::AtomTable) = at

@inline function Tables.getcolumn(at::AtomTable, nm::Symbol)
    col = Tables.getcolumn(_atoms(at), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _atoms(at)._idx_map[idx], getfield(at, :_idx))
    )
end

@inline function Base.getproperty(at::AtomTable, nm::Symbol)
    hasfield(typeof(at), nm) && return getfield(at, nm)
    Tables.getcolumn(at, nm)
end

@inline Tables.getcolumn(at::AtomTable, i::Int) = Tables.getcolumn(at, Tables.columnnames(at)[i])
@inline Tables.columnnames(at::AtomTable) = Tables.columnnames(_atoms(at))
@inline Tables.schema(at::AtomTable) = Tables.schema(_atoms(at))

@inline Base.size(at::AtomTable) = (length(getfield(at, :_idx)), length(_atom_table_cols))
@inline Base.size(at::AtomTable, dim) = size(at)[dim]
@inline Base.length(at::AtomTable) = size(at, 1)

function Base.push!(at::AtomTable{T}, t::AtomTuple{T}; kwargs...) where T
    sys = getfield(at, :_sys)
    push!(sys, t; kwargs...)
    push!(getfield(at, :_idx), sys._curr_idx)
    at
end

@inline function _filter_atoms(f::Function, sys::System{T}) where T
    AtomTable(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._atoms),
        :idx
    )))
end

@inline function Base.filter(f::Function, at::AtomTable)
    AtomTable(getfield(at, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, at),
        :idx
    )))
end

@inline function Base.iterate(at::AtomTable, st = 1)
    st > length(at) ?
        nothing :
        (atom_by_idx(getfield(at, :_sys), getfield(at, :_idx)[st]), st + 1)
end
@inline Base.eltype(::AtomTable{T}) where T = Atom{T}
@inline Base.getindex(at::AtomTable{T}, i::Int) where T = atom_by_idx(getfield(at, :_sys), getfield(at, :_idx)[i])
@inline Base.keys(at::AtomTable) = LinearIndices((length(at),))

"""
    $(TYPEDEF)

Mutable representation of an individual atom in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `element::ElementType`
 - `name::String`
 - `atom_type::String`
 - `r::Vector3{T}`
 - `v::Vector3{T}`
 - `F::Vector3{T}`
 - `formal_charge::Int`
 - `charge::T`
 - `radius::T`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Atom(
    number::Int,
    element::ElementType,
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags();
    # keyword arguments
    frame_id::Int = 1
)
```
Creates a new `Atom{Float32}` in the default system.

```julia
Atom(
    ac::AbstractAtomContainer{T},
    number::Int,
    element::ElementType,
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags();
    # keyword arguments
    frame_id::Int = 1
)
```
Creates a new `Atom{T}` in the given atom container (e.g., `System{T}` or `Molecule{T}`).

    Atom(a::AtomTuple{T}; frame_id::Int = 1)
    Atom(sys::System{T}, a::AtomTuple{T}; frame_id::Int = 1)

Constructor variants that create a new system atom based on the given `AtomTuple{T}`. The new atom
is automatically assigned a new `idx`.
"""
@auto_hash_equals struct Atom{T} <: AbstractSystemComponent{T}
    _sys::System{T}
    _row::_AtomTableRow{T}
end

function Atom(
    sys::System{T},
    number::Int,
    element::ElementType,
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags();
    frame_id::Int = 1,
    # private kwargs
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing,
    fragment_id::MaybeInt = nothing,
    nucleotide_id::MaybeInt = nothing,
    residue_id::MaybeInt = nothing
) where T
    idx = _next_idx(sys)
    push!(sys._atoms, AtomTuple{T}(number, element;
            idx = idx,
            name = name,
            atom_type = atom_type,
            r = r,
            v = v,
            F = F,
            formal_charge = formal_charge,
            charge = charge,
            radius = radius,
            properties = properties,
            flags = flags
        );
        frame_id = frame_id,
        molecule_id = molecule_id,
        chain_id = chain_id,
        fragment_id = fragment_id,
        nucleotide_id = nucleotide_id,
        residue_id = residue_id
    )
    atom_by_idx(sys, idx)::Atom{T}
end

function Atom(
    ac::AbstractAtomContainer{T},
    number::Int,
    element::ElementType,
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags();
    frame_id::Int = 1
) where T
    push!(ac, 
         (idx=0, number=number, element=element, name=name, 
          atom_type=atom_type, r=r, v=v, F=F, formal_charge=formal_charge, 
          charge=charge, radius=radius, properties=properties, flags=flags)::AtomTuple{T};
          frame_id=frame_id)
    atom_by_idx(ac._sys, ac._sys._curr_idx)::Atom{T}
end

function Atom(
    number::Int,
    element::ElementType,
    name::String = "",
    atom_type::String = "",
    r::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    v::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    F::Vector3{Float32} = Vector3{Float32}(0, 0, 0),
    formal_charge::Int = 0,
    charge::Float32 = 0.0f32,
    radius::Float32 = 0.0f32,
    properties::Properties = Properties(),
    flags::Flags = Flags();
    kwargs...
)
    Atom(default_system(), number, element, name, atom_type, r, v, F, formal_charge, charge, radius,
        properties, flags; kwargs...)::Atom{Float32}
end

@inline function Atom(sys::System{T}, atom::AtomTuple{T}; kwargs...) where T
    Atom(sys, ntuple(i -> atom[i+1], length(atom) - 1)...; kwargs...)::Atom{T}
end

@inline function Atom(t::AtomTuple{Float32}; kwargs...)
    Atom(default_system(), t; kwargs...)::Atom{Float32}
end

@inline Tables.rows(at::AtomTable) = at
@inline Tables.getcolumn(atom::Atom, name::Symbol) = Tables.getcolumn(getfield(atom, :_row), name)

@inline function Base.getproperty(atom::Atom, name::Symbol)
    hasfield(typeof(atom), name) && return getfield(atom, name)
    getproperty(getfield(atom, :_row), name)
end

@inline function Base.setproperty!(atom::Atom, name::Symbol, val)
    hasfield(typeof(atom), name) && return setfield!(atom, name, val)
    setproperty!(getfield(atom, :_row), name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", atom::Atom) = show(io, getfield(atom, :_row))
@inline Base.show(io::IO, atom::Atom) = show(io, getfield(atom, :_row))

@inline Base.parent(atom::Atom) = atom._sys
@inline parent_system(atom::Atom) = parent(atom)

@inline function parent_molecule(atom::Atom) 
    isnothing(atom._row.molecule_id) ?
        nothing :
        molecule_by_idx(parent(atom), atom._row.molecule_id)
end

@inline function parent_chain(atom::Atom)
    isnothing(atom._row.chain_id) ?
        nothing :
        chain_by_idx(atom._sys, atom._row.chain_id)
end

@inline function parent_fragment(atom::Atom)
    isnothing(atom._row.fragment_id) ?
        nothing :
        fragment_by_idx(parent(atom), atom._row.fragment_id)
end

@inline function parent_nucleotide(atom::Atom)
    isnothing(atom._row.nucleotide_id) ?
        nothing :
        nucleotide_by_idx(parent(atom), atom._row.nucleotide_id)
end

@inline function parent_residue(atom::Atom)
    isnothing(atom._row.residue_id) ?
        nothing :
        residue_by_idx(parent(atom), atom._row.residue_id)
end

"""
    $(TYPEDSIGNATURES)

Returns the `Atom{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
atom exists.
"""
@inline function atom_by_idx(sys::System{T}, idx::Int) where T
    Atom{T}(sys, _row_by_idx(sys._atoms, idx))
end

"""
    $(TYPEDSIGNATURES)

Returns the first `Atom{T}` associated with the given `name` in `ac`. Returns nothing if no such
atom exists.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.
"""
@inline function atom_by_name(
    ac::AbstractAtomContainer{T},
    name::String;
    frame_id::MaybeInt = 1
) where T
    idx = filter(atom -> atom.name == name, atoms(ac; frame_id = frame_id)).idx
    isempty(idx) ? nothing : atom_by_idx(parent(ac), first(idx))
end

"""
    atoms(::Chain)
    atoms(::Fragment)
    atoms(::Molecule)
    atoms(::Nucleotide)
    atoms(::Protein)
    atoms(::Residue)
    atoms(::System)

Returns an `AtomTable{T}` containing all atoms of the given atom container.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.
 - `molecule_id::Union{MaybeInt, Some{Nothing}} = nothing`: \
Any value other than `nothing` limits the results to atoms belonging to the given molecule ID.
 - `chain_id::Union{MaybeInt, Some{Nothing}} = nothing`: \
Any value other than `nothing` limits the results to atoms belonging to the given chain ID.
- `fragment_id::Union{MaybeInt, Some{Nothing}} = nothing`: \
Any value other than `nothing` limits the results to atoms belonging to the given fragment ID.
- `nucleotide_id::Union{MaybeInt, Some{Nothing}} = nothing`: \
Any value other than `nothing` limits the results to atoms belonging to the given nucleotide ID.
- `residue_id::Union{MaybeInt, Some{Nothing}} = nothing`: \
Any value other than `nothing` limits the results to atoms belonging to the given residue ID.
"""
@inline function atoms(sys::System{T};
    frame_id::MaybeInt = 1,
    molecule_id::Union{MaybeInt, Some{Nothing}} = nothing,
    chain_id::Union{MaybeInt, Some{Nothing}} = nothing,
    fragment_id::Union{MaybeInt, Some{Nothing}} = nothing,
    nucleotide_id::Union{MaybeInt, Some{Nothing}} = nothing,
    residue_id::Union{MaybeInt, Some{Nothing}} = nothing
) where T
    _filter_atoms(atom ->
        (isnothing(frame_id)      || atom.frame_id == frame_id) &&
        (isnothing(molecule_id)   || atom.molecule_id == something(molecule_id)) &&
        (isnothing(chain_id)      || atom.chain_id == something(chain_id)) &&
        (isnothing(fragment_id)   || atom.fragment_id == something(fragment_id)) &&
        (isnothing(nucleotide_id) || atom.nucleotide_id == something(nucleotide_id)) &&
        (isnothing(residue_id)    || atom.residue_id == something(residue_id)),
        sys
    )
end

"""
    atoms_df(::Chain)
    atoms_df(::Fragment)
    atoms_df(::Molecule)
    atoms_df(::Nucleotide)
    atoms_df(::Protein)
    atoms_df(::Residue)
    atoms_df(::System)

Returns a `DataFrame` containing all atoms of the given atom container.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.
"""
@inline function atoms_df(sys::System{T}; kwargs...) where T
    DataFrame(atoms(sys; kwargs...))
end

"""
    eachatom(::Chain)
    eachatom(::Fragment)
    eachatom(::Molecule)
    eachatom(::Nucleotide)
    eachatom(::Protein)
    eachatom(::Residue)
    eachatom(::System)

Returns an `Atom{T}` generator for all atoms of the given atom container.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.

# Example
```julia
for atom in eachatom(sys)
    println(atom.name)
end
```
"""
@inline function eachatom(sys::System{T}; kwargs...) where T
    (atom for atom in atoms(sys; kwargs...))
end

"""
    natoms(::Chain)
    natoms(::Fragment)
    natoms(::Molecule)
    natoms(::Nucleotide)
    natoms(::Protein)
    natoms(::Residue)
    natoms(::System)

Returns the number of atoms in the given atom container.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.
"""
@inline function natoms(sys::System; kwargs...)
    length(atoms(sys; kwargs...))
end

"""
    bonds(::Atom)

Returns a `BondTable{T}` containing all bonds of the given atom.
"""
@inline function bonds(atom::Atom)
    _filter_bonds(
        bond -> bond.a1 == atom.idx || bond.a2 == atom.idx,
        parent(atom)
    )
end

"""
    bonds_df(::Atom)

Returns a `SubDataFrame` containing all bonds of the given atom.
"""
@inline function bonds_df(atom::Atom)
    DataFrame(bonds(atom))
end

"""
    eachbond(::Atom)

Returns a `Bond{T}` generator for all bonds of the given atom.
"""
@inline function eachbond(atom::Atom{T}) where T
    (bond for bond in bonds(atom))
end

"""
    nbonds(::Atom)

Returns the number of bonds of the given atom.
"""
@inline function nbonds(atom::Atom)
    length(bonds(atom))
end

"""
    push!(::Fragment{T}, atom::AtomTuple{T})
    push!(::Molecule{T}, atom::AtomTuple{T})
    push!(::Nucleotide{T}, atom::AtomTuple{T})
    push!(::Protein{T}, atom::AtomTuple{T})
    push!(::Residue{T}, atom::AtomTuple{T})
    push!(::System{T}, atom::AtomTuple{T})

Creates a new atom in the given atom container, based on the given tuple. The new atom is
automatically assigned a new `idx`.

# Supported keyword arguments
 - `frame_id::Int = 1`
"""
@inline function Base.push!(sys::System{T}, atom::AtomTuple{T};
    frame_id::Int = 1,
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing,
    fragment_id::MaybeInt = nothing,
    nucleotide_id::MaybeInt = nothing,
    residue_id::MaybeInt = nothing
) where T
    push!(sys._atoms,
        (; atom..., idx = _next_idx(sys));
        frame_id = frame_id,
        molecule_id = molecule_id,
        chain_id = chain_id,
        fragment_id = fragment_id,
        nucleotide_id = nucleotide_id,
        residue_id = residue_id
    )
    sys
end

@enumx FullNameType begin
    # Do not add extensions
    NO_VARIANT_EXTENSIONS = 1
    # Add the residue extensions
    ADD_VARIANT_EXTENSIONS = 2
    # Add the residue ID
    ADD_RESIDUE_ID = 3
    # Add the residue ID and the residue extension
    ADD_VARIANT_EXTENSIONS_AND_ID = 4
end

function get_full_name(
    a::Atom{T},
    type::FullNameType.T = FullNameType.ADD_VARIANT_EXTENSIONS
) where {T<:Real}

    # determine the parent's name
    f = parent_fragment(a)

    parent_name = ""

    if isnothing(f)
        # look for a molecule containing the atom
        m = parent_molecule(a)
        parent_name = strip(m.name)
    else
        # retrieve the fragment name
        parent_name = get_full_name(f, type)
    end

    # retrieve the atom name
    name = strip(a.name)

    # add the parent name only if non-empty
    if !isempty(parent_name)
        name = string(parent_name, ":", name)
    end

    name
end

@inline distance(a1::Atom, a2::Atom) = distance(a1.r, a2.r)

"""
    $(TYPEDSIGNATURES)

    Decides if two atoms are bound to each other.
    Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.
"""
function is_bound_to(a1::Atom, a2::Atom)
    s = parent(a1)

    if s != parent(a2)
        return false
    end

    return !isnothing(
        findfirst(
            b::Bond -> 
                ((b.a1 == a1.idx) && (b.a2 == a2.idx)) ||
                ((b.a1 == a2.idx) && (b.a2 == a1.idx)), 
            non_hydrogen_bonds(s)
        )
    )
end

"""
    $(TYPEDSIGNATURES)

    Decides if two atoms are geminal.
    
    Two atoms are geminal if they do not share a common bond but both have a
    bond to a third atom. For example the two hydrogen atoms in water are geminal. 
    Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.
"""
function is_geminal(a1::Atom, a2::Atom)
    if a1 == a2
        return false
    end

    # an atom is geminal to another, if it is not directly bonded to it...
    is_geminal = !is_bound_to(a1, a2)

    # ...and is bonded to an atom that is bonded to the other atom
    is_geminal && any(map(b -> is_bound_to(get_partner(b, a1), a2), non_hydrogen_bonds(a1)))
end

"""
    $(TYPEDSIGNATURES)

Decides if two atoms are vicinal.

Two atoms are vicinal if they are separated by three bonds (1-4 position).
Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.
"""
function is_vicinal(a1::Atom, a2::Atom)
    if a1 == a2
        return false
    end

    # an atom is vicinal to another, if it is not directly bonded to it...
    is_vicinal = !is_bound_to(a1, a2)

    # ...and is bonded to an atom that is bonded to an atom that is bonded to this atom
    if is_vicinal
        is_vicinal = false

        for b_1 in non_hydrogen_bonds(a1)
            partner_1 = get_partner(b_1, a1)

            for b_2 in non_hydrogen_bonds(partner_1)
                partner_2 = get_partner(b_2, partner_1)

                if is_bound_to(partner_2, a2)
                    return true
                end
            end
        end
    end

    return false
end
