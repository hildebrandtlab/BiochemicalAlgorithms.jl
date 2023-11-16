export Substructure, filter_atoms

@auto_hash_equals struct Substructure{T<:Real} <: AbstractMolecule{T}
    name::String

    parent::AbstractAtomContainer{T}
    
    _atoms::SubDataFrame
    _bonds::SubDataFrame
    
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

function filter_atoms(fn, mol; name="", adjacent_bonds=false)
    atom_view = filter(fn, _atoms(mol), view=true)
    bond_view = filter(
        [:a1, :a2] => (a1, a2) -> 
            adjacent_bonds ? a1 ∈ atom_view.idx || a2 ∈ atom_view.idx
                           : a1 ∈ atom_view.idx && a2 ∈ atom_view.idx,
        _bonds(mol), view=true)

    Substructure(name, mol, atom_view, bond_view)
end

function Base.copy(substruct::Substructure{T}) where T
    sys = System{T}(substruct.name)
    sys._curr_idx = sys._curr_idx

    sys.properties = copy(substruct.properties)
    sys.flags      = copy(substruct.parent.flags)

    sys._atoms = IndexedDataFrame(copy(substruct._atoms))
    sys._bonds = copy(substruct._bonds)

    sys._molecules   = copy(_molecules(substruct))
    sys._chains      = copy(_chains(substruct))
    sys._fragments   = copy(_fragments(substruct))
    sys._nucleotides = copy(_nucleotides(substruct))
    sys._residues    = copy(_residues(substruct))

    sys
end

"""
$(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's atoms matching the given criteria (value or
`missing`). Fields given as `nothing` are ignored. The returned `DataFrame` contains all public and
private atom fields.
"""
function _atoms(substruct::Substructure{T};
    frame_id::Union{Nothing, Int} = 1,
    molecule_id::Union{Nothing, MaybeInt} = nothing,
    chain_id::Union{Nothing, MaybeInt} = nothing,
    fragment_id::Union{Nothing, MaybeInt} = nothing,
    nucleotide_id::Union{Nothing, MaybeInt} = nothing,
    residue_id::Union{Nothing, MaybeInt} = nothing
) where T
    cols = Tuple{Symbol, MaybeInt}[]
    isnothing(frame_id)      || push!(cols, (:frame_id, frame_id))
    isnothing(molecule_id)   || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)      || push!(cols, (:chain_id, chain_id))
    isnothing(fragment_id)   || push!(cols, (:fragment_id, fragment_id))
    isnothing(nucleotide_id) || push!(cols, (:nucleotide_id, nucleotide_id))
    isnothing(residue_id)    || push!(cols, (:residue_id, residue_id))

    get(
        groupby(substruct._atoms, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        view(substruct._atoms, Int[], :)
    )::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int}}
end

function _bonds(substruct::Substructure; kwargs...)
    aidx = _atoms(substruct; kwargs...).idx
    @rsubset(
        substruct._bonds, :a1 in aidx || :a2 in aidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _molecules(substruct::Substructure; kwargs...)
    midx = unique(_atoms(substruct; kwargs...).molecule_id)
    @rsubset(
        substruct.parent._molecules, :idx in midx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _chains(substruct::Substructure; kwargs...)
    cidx = _atoms(substruct; kwargs...).chain_id
    @rsubset(
        substruct.parent._chains, :idx in cidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _fragments(substruct::Substructure; kwargs...)
    fidx = _atoms(substruct; kwargs...).fragment_id
    @rsubset(
        substruct.parent._fragments, :idx in fidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _nucleotides(substruct::Substructure; kwargs...)
    nidx = _atoms(substruct; kwargs...).nucleotide_id
    @rsubset(
        substruct.parent._nucleotides, :idx in nidx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

function _residues(substruct::Substructure; kwargs...)
    ridx = _atoms(substruct; kwargs...).residue_id
    @rsubset(
        substruct.parent._residues, :idx in ridx; view = true
    )::SubDataFrame{DataFrame, DataFrames.Index, <:AbstractVector{Int}}
end

@inline function eachatom(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Atom{T}(sys, row) for row in eachrow(_atoms(substruct; kwargs...)))
end

@inline function atoms(substruct::Substructure; kwargs...)
    collect(eachatom(substruct; kwargs...))
end

@inline function eachbond(substruct::Substructure{T}; kwargs...) where T
    sys = substruct.parent isa System{T} ? substruct.parent : parent_system(substruct.parent)
    (Bond{T}(sys, row) for row in eachrow(_bonds(substruct; kwargs...)))
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
    view(_atoms(ac; kwargs...), :, 1:length(fieldnames(AtomTuple{T})))
end

function bonds_df(ac::Substructure{T}; kwargs...) where {T<:Real}
    _bonds(ac; kwargs...)
end

@inline function natoms(substruct::Substructure; kwargs...)
    nrow(_atoms(substruct; kwargs...))
end

@inline function nbonds(substruct::Substructure; kwargs...)
    nrow(_bonds(substruct; kwargs...))
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
