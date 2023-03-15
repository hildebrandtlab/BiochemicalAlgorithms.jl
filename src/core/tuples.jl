export AtomTuple, BondTuple, MoleculeTuple, ChainTuple, FragmentTuple, NucleotideTuple, ResidueTuple

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

@inline function _with_idx(atom::AtomTuple{T}, idx::Int) where T
    ntuple(i -> i == 1 ? idx : atom[i], length(atom))
end

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

const MoleculeTuple = @NamedTuple begin
    idx::Int
    name::String
    properties::Properties
end

@inline function _with_idx(mol::MoleculeTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : mol[i], length(mol))
end

const ChainTuple = MoleculeTuple

const FragmentTuple = @NamedTuple begin
    idx::Int
    number::Int
    name::String
    properties::Properties
end

@inline function _with_idx(frag::FragmentTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end

const NucleotideTuple = FragmentTuple

const ResidueTuple = @NamedTuple begin
    idx::Int
    number::Int
    type::AminoAcid
    properties::Properties
end

@inline function _with_idx(frag::ResidueTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end
