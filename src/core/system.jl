export System

### Types
struct System{T<:Real}
    molecules::Vector{Molecule{T}}
end

System() = System{Float32}()

### Functions
function Base.push!(s::System, m::Molecule)
    push!(s.molecules, m)
end
