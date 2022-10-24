export Atom

const Atom{T} = @NamedTuple begin
    number::Int
    name::String
    element::Element
    atomtype::String
    r::Vector3{T}
    v::Vector3{T}
    F::Vector3{T}
    has_velocity::Bool
    has_force::Bool
    frame_id::Int
    props::Properties
end
 
Atom() = Atom{Float32}()
