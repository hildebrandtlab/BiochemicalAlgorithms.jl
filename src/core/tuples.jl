export AtomTuple, BondTuple, MoleculeTuple, ChainTuple, FragmentTuple, NucleotideTuple, ResidueTuple

"""
    const AtomTuple{T} = NamedTuple{...}

Tuple-based atom representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `element::ElementType`
 - `name::String`
 - `atomtype::String`
 - `r::Vector3{T}`
 - `v::Vector3{T}`
 - `F::Vector3{T}`
 - `has_velocity::Bool`
 - `has_force::Bool`
 - `properties::Properties`
"""
const AtomTuple{T} = @NamedTuple begin
    idx::Int
    number::Int
    element::ElementType
    name::String
    atomtype::String
    r::Vector3{T}
    v::Vector3{T}
    F::Vector3{T}
    has_velocity::Bool
    has_force::Bool
    properties::Properties
end

"""
    _with_idx(::AtomTuple{T}, idx::Int)
    _with_idx(::BondTuple{T}, idx::Int)
    _with_idx(::ChainTuple{T}, idx::Int)
    _with_idx(::FragmentTuple{T}, idx::Int)
    _with_idx(::MoleculeTuple{T}, idx::Int)
    _with_idx(::NucleotideTuple{T}, idx::Int)
    _with_idx(::ProteinTuple{T}, idx::Int)
    _with_idx(::ResidueTuple{T}, idx::Int)

Returns a copy of the given tuple with replaced `idx`.
"""
@inline function _with_idx(atom::AtomTuple{T}, idx::Int) where T
    ntuple(i -> i == 1 ? idx : atom[i], length(atom))
end

"""
    const BondTuple = NamedTuple{...}

Tuple-based bond representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `a1::Int`
 - `a2::Int`
 - `order::BondOrderType`
 - `properties::Properties`
"""
const BondTuple = @NamedTuple begin
    idx::Int
    a1::Int
    a2::Int
    order::BondOrderType
    properties::Properties
end

@inline function _with_idx(bond::BondTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : bond[i], length(bond))
end

"""
    const MoleculeTuple = NamedTuple{...}

Tuple-based molecule representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`
"""
const MoleculeTuple = @NamedTuple begin
    idx::Int
    name::String
    properties::Properties
end

@inline function _with_idx(mol::MoleculeTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : mol[i], length(mol))
end

"""
    const ChainTuple = NamedTuple{...}

Tuple-based chain representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`
"""
const ChainTuple = MoleculeTuple

"""
    const FragmentTuple = NamedTuple{...}

Tuple-based fragment representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`
"""
const FragmentTuple = @NamedTuple begin
    idx::Int
    number::Int
    name::String
    properties::Properties
end

@inline function _with_idx(frag::FragmentTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end

"""
    const NucleotideTuple = NamedTuple{...}

Tuple-based nucleotide representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`
"""
const NucleotideTuple = FragmentTuple

"""
    const ResidueTuple = NamedTuple{...}

Tuple-based residue representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `type::AminoAcid`
 - `properties::Properties`
"""
const ResidueTuple = @NamedTuple begin
    idx::Int
    number::Int
    type::AminoAcid
    properties::Properties
end

@inline function _with_idx(frag::ResidueTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end
