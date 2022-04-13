export Atom

const Atom = @NamedTuple begin
    id::Int
    local_id::Int
    element::Element
    atomtype::String
end