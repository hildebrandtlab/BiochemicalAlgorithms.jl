export Atom, AtomProperties

const Atom{T} = @NamedTuple begin
    number::Int
    name::String
    element::ElementType
    atomtype::String
    r::Vector3{T}
    v::Vector3{T}
    F::Vector3{T}
    has_velocity::Bool
    has_force::Bool
    frame_id::Int
end

const AtomProperties = Dict{String, Any}
