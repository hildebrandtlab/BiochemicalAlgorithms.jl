export
    Substructure,
    filter_atoms

@auto_hash_equals struct Substructure{T<:Real, A<:AbstractAtomContainer{T}} <: AbstractAtomContainer{T}
    name::String

    parent::A

    _atoms::AtomTable{T}
    _bonds::BondTable{T}

    properties::Properties

    function Substructure{T,A}(
        name::AbstractString,
        parent::A,
        atoms::AtomTable{T},
        bonds::BondTable{T},
        properties::Properties = parent.properties
    ) where {T<:Real, A<:AbstractAtomContainer{T}}
        new(name, parent, atoms, bonds, properties)
    end
end

@inline Substructure(
    name::AbstractString,
    parent::A,
    atoms::AtomTable{T},
    bonds::BondTable{T},
    properties::Properties = parent.properties
) where {T, A} = Substructure{T, typeof(parent)}(name, parent, atoms, bonds, properties)

function filter_atoms(fn, mol::AbstractAtomContainer; name::AbstractString="", adjacent_bonds::Bool=false)
    atom_view = filter(fn, atoms(mol))
    idxset = Set(atom_view.idx)
    bond_view = filter(b ->
        adjacent_bonds ? b.atom1_idx ∈ idxset || b.atom2_idx ∈ idxset
                       : b.atom1_idx ∈ idxset && b.atom2_idx ∈ idxset,
        bonds(parent(mol))
    )
    Substructure(name, mol, atom_view, bond_view)
end

function Base.copy(substruct::Substructure{T}) where T
    sys = System{T}(substruct.name)
    sys._curr_idx = parent(substruct)._curr_idx

    sys.properties = copy(substruct.properties)
    sys.flags      = copy(substruct.parent.flags)

    sys._atoms = _atom_table(T, deepcopy(substruct._atoms))
    sys._bonds = _bond_table(deepcopy(substruct._bonds))

    sys._molecules   = _molecule_table(deepcopy(molecules(substruct)))
    sys._chains      = _chain_table(deepcopy(chains(substruct)))
    sys._secondary_structures = _secondary_structure_table(deepcopy(secondary_structures(substruct)))
    sys._fragments   = _fragment_table(deepcopy(fragments(substruct)))

    sys
end

"""
$(TYPEDSIGNATURES)

Returns an `AtomTable` for all of the given system's atoms matching the given criteria (value or
`missing`). Fields given as `nothing` are ignored. The returned table contains all public and
private atom fields.
"""
@inline function atoms(substruct::Substructure{T};
    frame_id::MaybeInt = 1,
    molecule_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    chain_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    fragment_idx::Union{MaybeInt, Some{Nothing}} = nothing
) where T
    filter(row ->
        (isnothing(frame_id)                || row.frame_id == frame_id) &&
        (isnothing(molecule_idx)            || row.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)               || row.chain_idx == something(chain_idx)) &&
        (isnothing(fragment_idx)            || row.fragment_idx == something(fragment_idx)),
        substruct._atoms
    )
end

@inline function bonds(substruct::Substructure; kwargs...)
    aidx = Set(atoms(substruct; kwargs...).idx)
    filter(row -> row.atom1_idx in aidx || row.atom2_idx in aidx, substruct._bonds)
end

@inline function molecules(substruct::Substructure; kwargs...)
    midx = Set(atoms(substruct; kwargs...).molecule_idx)
    filter(row -> row.idx in midx, molecules(parent(substruct)))
end

@inline function proteins(substruct::Substructure; kwargs...)
    midx = Set(atoms(substruct; kwargs...).molecule_idx)
    filter(row -> row.idx in midx, proteins(parent(substruct)))
end

@inline function chains(substruct::Substructure; kwargs...)
    cidx = Set(atoms(substruct; kwargs...).chain_idx)
    filter(row -> row.idx in cidx, chains(parent(substruct)))
end

@inline function secondary_structures(substruct::Substructure{T}; kwargs...) where T
    fidx = setdiff(Set(atoms(substruct; kwargs...).fragment_idx), Set([nothing]))
    ft = FragmentTable{T}(parent(substruct), collect(Int, fidx))
    secondary_structures(ft)
end

@inline function fragments(substruct::Substructure; kwargs...)
    fidx = Set(atoms(substruct; kwargs...).fragment_idx)
    filter(row -> row.idx in fidx, fragments(parent(substruct)))
end

@inline function residues(substruct::Substructure; kwargs...)
    fidx = Set(atoms(substruct; kwargs...).fragment_idx)
    filter(row -> row.idx in fidx, residues(parent(substruct)))
end

@inline function nucleotides(substruct::Substructure; kwargs...)
    fidx = Set(atoms(substruct; kwargs...).fragment_idx)
    filter(row -> row.idx in fidx, nucleotides(parent(substruct)))
end

# natoms, nbonds, ...
for comp in (
    :atoms, :bonds, :molecules, :proteins, :chains,
    :secondary_structures, :fragments, :residues, :nucleotides
)
    @eval begin
        @inline function $(Symbol("n", comp))(substruct::Substructure; kwargs...)
            length($(comp)(substruct; kwargs...))
        end
    end
end

@inline Base.parent(substruct::Substructure) = parent(substruct.parent)
@inline parent_system(substruct::Substructure) = parent_system(substruct.parent)

@inline function atom_by_idx(substruct::Substructure{T}, idx::Int) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    atom_by_idx(sys, idx)
end
