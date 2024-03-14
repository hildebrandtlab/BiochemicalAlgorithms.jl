export Substructure, filter_atoms

@auto_hash_equals struct Substructure{T<:Real, A<:AbstractAtomContainer{T}} <: AbstractMolecule{T}
    name::String

    parent::A
    
    _atoms::AtomTable{T}
    _bonds::BondTable{T}
    
    properties::Properties

    function Substructure{T,A}(
        name::String,
        parent::A,
        atoms::AtomTable{T},
        bonds::BondTable{T},
        properties::Properties = parent.properties
    ) where {T<:Real, A<:AbstractAtomContainer{T}}
        new(name, parent, atoms, bonds, properties)
    end
end

@inline Substructure(
    name::String,
    parent::A,
    atoms::AtomTable{T},
    bonds::BondTable{T},
    properties::Properties = parent.properties
) where {T, A} = Substructure{Float32, typeof(parent)}(name, parent, atoms, bonds, properties)

function filter_atoms(fn, mol::AbstractAtomContainer{T}; name="", adjacent_bonds=false) where T
    atom_view = filter(fn, atoms(mol))
    idxset = Set(atom_view.idx)
    bond_view = filter(row ->
        adjacent_bonds ? row.a1 ∈ idxset || row.a2 ∈ idxset
                        : row.a1 ∈ idxset && row.a2 ∈ idxset,
        bonds(mol)
    )
    Substructure(name, mol, atom_view, bond_view)
end

function Base.copy(substruct::Substructure{T}) where T
    sys = System{T}(substruct.name)
    sys._curr_idx = sys._curr_idx

    sys.properties = copy(substruct.properties)
    sys.flags      = copy(substruct.parent.flags)

    sys._atoms = _atom_table(T, deepcopy(substruct._atoms))
    sys._bonds = _bond_table(deepcopy(substruct._bonds))

    sys._molecules   = _molecule_table(deepcopy(molecules(substruct)))
    sys._chains      = _chain_table(deepcopy(chains(substruct)))
    sys._fragments   = _fragment_table(deepcopy(fragments(substruct)))
    sys._nucleotides = _nucleotide_table(deepcopy(nucleotides(substruct)))
    sys._residues    = _residue_table(deepcopy(residues(substruct)))

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
    fragment_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    nucleotide_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    residue_idx::Union{MaybeInt, Some{Nothing}} = nothing
) where T
    filter(row ->
        (isnothing(frame_id)       || row.frame_id == frame_id) &&
        (isnothing(molecule_idx)   || row.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)      || row.chain_idx == something(chain_idx)) &&
        (isnothing(fragment_idx)   || row.fragment_idx == something(fragment_idx)) &&
        (isnothing(nucleotide_idx) || row.nucleotide_idx == something(nucleotide_idx)) &&
        (isnothing(residue_idx)    || row.residue_idx == something(residue_idx)),
        substruct._atoms
    )
end

@inline function bonds(substruct::Substructure; kwargs...)
    aidx = Set(atoms(substruct; kwargs...).idx)
    filter(row -> row.a1 in aidx || row.a2.idx, substruct._bonds)
end

@inline function molecules(substruct::Substructure; kwargs...)
    midx = Set(atoms(substruct; kwargs...).molecule_idx)
    filter(row -> row.idx in midx, molecules(substruct.parent))
end

@inline function chains(substruct::Substructure; kwargs...)
    cidx = Set(atoms(substruct; kwargs...).chain_idx)
    filter(row -> row.idx in cidx, chains(substruct.parent))
end

@inline function fragments(substruct::Substructure; kwargs...)
    fidx = Set(atoms(substruct; kwargs...).fragment_idx)
    filter(row -> row.idx in fidx, fragments(substruct.parent))
end

@inline function nucleotides(substruct::Substructure; kwargs...)
    nidx = Set(atoms(substruct; kwargs...).nucleotide_idx)
    filter(row -> row.idx in nidx, nucleotides(substruct.parent))
end

@inline function residues(substruct::Substructure; kwargs...)
    ridx = Set(atoms(substruct; kwargs...).residue_idx)
    filter(row -> row.idx in ridx, residues(substruct.parent))
end

@inline function natoms(substruct::Substructure; kwargs...)
    length(atoms(substruct; kwargs...))
end

@inline function nbonds(substruct::Substructure; kwargs...)
    length(bonds(substruct; kwargs...))
end

@inline function nfragments(substruct::Substructure; kwargs...)
    length(fragments(substruct; kwargs...))
end

@inline Base.parent(substruct::Substructure) = parent(substruct.parent)
@inline parent_system(substruct::Substructure) = parent_system(substruct.parent)

@inline function atom_by_idx(substruct::Substructure{T}, idx::Int) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    atom_by_idx(sys, idx)
end
