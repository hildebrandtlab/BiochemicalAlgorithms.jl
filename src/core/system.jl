using DataFrames

export System, add_molecule

### Types
struct System{T<:Real}
    molecules::DataFrame
    atoms::DataFrame
    positions::DataFrame
    velocities::DataFrame
    forces::DataFrame

    function System{T}() where {T<:Real}
        molecules = DataFrame(id=Int[], name=String[])
        
        atoms = DataFrame(Atom[])

        positions = DataFrame(id=Int[], atom_id=Int[], r=Vector3{T}[])

        velocities = DataFrame(id=Int[], atom_id=Int[], v=Vector3{T}[])

        forces = DataFrame(id=Int[], atom_id=Int[], f=Vector3{T}[])

        new(molecules, atoms, positions, velocities, forces)
    end
end

System() = System{Float32}()

### Functions

function add_molecule(s::System, name::String)
    new_id = isempty(s.molecules) ? 1 : maximum(s.molecules.id)+1

    mol = Molecule(new_id, name)
    push!(s.molecules, [mol.id, mol.name])

    return mol
end