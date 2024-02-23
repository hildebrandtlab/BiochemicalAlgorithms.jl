export Substructure, filter_atoms

@auto_hash_equals struct Substructure{T<:Real} <: AbstractMolecule{T}
    name::String

    parent::AbstractAtomContainer{T}
    
    _atoms::_AtomTable{T}
    _bonds::_BondTable
    
    properties::Properties

    function Substructure{T}(name,
                             parent, 
                             atoms, 
                             bonds, 
                             properties = parent.properties) where {T<:Real}
        new(name, parent, atoms, bonds, properties)
    end
end

Substructure(name,
             parent,
             atoms,
             bonds,
             properties = parent.properties) = 
                Substructure{Float32}(name, parent, atoms, bonds, properties)

function filter_atoms(fn, mol::AbstractAtomContainer{T}; name="", adjacent_bonds=false) where T
    atom_view = Tables.materializer(_AtomTable{T})(TableOperations.filter(fn, _atoms(mol)))
    idxset = Set(atom_view.idx)
    bond_view = Tables.materializer(_BondTable)(
        TableOperations.filter(row ->
            adjacent_bonds ? row.a1 ∈ idxset || row.a2 ∈ idxset
                           : row.a1 ∈ idxset && row.a2 ∈ idxset,
            _bonds(mol)
        )
    )
    Substructure(name, mol, atom_view, bond_view)
end

function Base.copy(substruct::Substructure{T}) where T
    sys = System{T}(substruct.name)
    sys._curr_idx = sys._curr_idx

    sys.properties = copy(substruct.properties)
    sys.flags      = copy(substruct.parent.flags)

    sys._atoms = deepcopy(substruct._atoms)
    sys._bonds = deepcopy(substruct._bonds)

    sys._molecules   = IndexedDataFrame(copy(_molecules(substruct)))
    sys._chains      = IndexedDataFrame(copy(_chains(substruct)))
    sys._fragments   = IndexedDataFrame(copy(_fragments(substruct)))
    sys._nucleotides = IndexedDataFrame(copy(_nucleotides(substruct)))
    sys._residues    = IndexedDataFrame(copy(_residues(substruct)))

    sys
end

"""
$(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's atoms matching the given criteria (value or
`missing`). Fields given as `nothing` are ignored. The returned `DataFrame` contains all public and
private atom fields.
"""
function _atoms(substruct::Substructure{T};
    frame_id::MaybeInt = 1,
    molecule_id::Union{MaybeInt, Some{Nothing}} = nothing,
    chain_id::Union{MaybeInt, Some{Nothing}} = nothing,
    fragment_id::Union{MaybeInt, Some{Nothing}} = nothing,
    nucleotide_id::Union{MaybeInt, Some{Nothing}} = nothing,
    residue_id::Union{MaybeInt, Some{Nothing}} = nothing
) where T
    TableOperations.filter(row ->
        (isnothing(frame_id)      || row.frame_id == frame_id) &&
        (isnothing(molecule_id)   || row.molecule_id == something(molecule_id)) &&
        (isnothing(chain_id)      || row.chain_id == something(chain_id)) &&
        (isnothing(fragment_id)   || row.fragment_id == something(fragment_id)) &&
        (isnothing(nucleotide_id) || row.nucleotide_id == something(nucleotide_id)) &&
        (isnothing(residue_id)    || row.residue_id == something(residue_id)),
        substruct._atoms
    )
end

function _bonds(substruct::Substructure; kwargs...)
    aidx = Set(_filter_select(_atoms(substruct; kwargs...), :idx))
    TableOperations.filter(row -> row.a1 in aidx || row.a2.idx, substruct._bonds)
end

function _molecules(substruct::Substructure; kwargs...)
    midx = Set(_filter_select(_atoms(substruct; kwargs...), :molecule_id))
    @rsubset(
        _molecules(substruct.parent), :idx in midx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _chains(substruct::Substructure; kwargs...)
    cidx = Set(_filter_select(_atoms(substruct; kwargs...), :chain_id))
    @rsubset(
        _chains(substruct.parent), :idx in cidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _fragments(substruct::Substructure; kwargs...)
    fidx = Set(_filter_select(_atoms(substruct; kwargs...), :fragment_id))
    @rsubset(
        _fragments(substruct.parent), :idx in fidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _nucleotides(substruct::Substructure; kwargs...)
    nidx = Set(_filter_select(_atoms(substruct; kwargs...), :nucleotide_id))
    @rsubset(
        _nucleotides(substruct.parent), :idx in nidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _residues(substruct::Substructure; kwargs...)
    ridx = Set(_filter_select(_atoms(substruct; kwargs...), :residue_id))
    @rsubset(
        _residues(substruct.parent), :idx in ridx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

@inline function eachatom(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Atom{T}(sys, row) for row in _atoms(substruct; kwargs...))
end

@inline function atoms(substruct::Substructure; kwargs...)
    collect(eachatom(substruct; kwargs...))
end

@inline function eachbond(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Bond{T}(sys, row) for row in _bonds(substruct; kwargs...))
end

@inline function bonds(substruct::Substructure; kwargs...)
    collect(eachbond(substruct; kwargs...))
end

@inline function eachfragment(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Fragment{T}(sys, row) for row in eachrow(_fragments(substruct; kwargs...)))
end

@inline function fragments(substruct::Substructure{T}; kwargs...) where T
    collect(eachfragment(substruct; kwargs...))
end

@inline function eachchain(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Chain{T}(sys, row) for row in eachrow(_chains(substruct; kwargs...)))
end

@inline function chains(substruct::Substructure{T}; kwargs...) where T
    collect(eachchain(substruct; kwargs...))
end

function atoms_df(ac::Substructure{T}; kwargs...) where {T<:Real}
    DataFrame(_atoms(ac; kwargs...))
end

function bonds_df(ac::Substructure{T}; kwargs...) where {T<:Real}
    DataFrame(_bonds(ac; kwargs...))
end

@inline function natoms(substruct::Substructure; kwargs...)
    count(_ -> true, _atoms(substruct; kwargs...))
end

@inline function nbonds(substruct::Substructure; kwargs...)
    count(_ -> true, _bonds(substruct; kwargs...))
end

@inline function nfragments(substruct::Substructure; kwargs...)
    nrow(_fragments(substruct; kwargs...))
end

@inline Base.parent(substruct::Substructure) = parent(substruct.parent)
@inline parent_system(substruct::Substructure) = parent_system(substruct.parent)

@inline function atom_by_idx(substruct::Substructure{T}, idx::Int) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    atom_by_idx(sys, idx)
end
