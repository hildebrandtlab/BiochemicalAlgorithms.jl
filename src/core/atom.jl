export Atom, has_property, get_property, set_property

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
    properties::Properties
end

function has_property(a::Atom, key::String)
    haskey(a.properties, key)
end

function get_property(a::Atom, key::String)
    a.properties[key]
end

function set_property(a::Atom, key::String, value)
    a.properties[key] = value
end
