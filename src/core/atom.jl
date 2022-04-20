export Atom

const Atom{T} = @NamedTuple begin
    id::Int
    molecule_id::Int
    frame_id::Int
    number::Int
    element::Element
    atomtype::String
    r::Vector3{T}
    v::Vector3{T}
    F::Vector3{T}
    has_velocity::Bool
    has_force::Bool
end