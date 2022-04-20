using DataFrames

export System, add_molecule!, add_atom!

### Types
struct System{T<:Real}
    molecules::DataFrame
    atoms::DataFrame

    function System{T}() where {T<:Real}
        molecules = DataFrame(id=Int[], name=String[])
        
        atoms = DataFrame(Atom{T}[])

        new(molecules, atoms)
    end
end

System() = System{Float32}()

### Functions

#TODO: for all functions:
#   - do we really need the id-column or is there a way to access row numbers?
#   - are these functions thread safe?
function add_molecule!(s::System, name::String)
    new_id = isempty(s.molecules) ? 1 : maximum(s.molecules.id)+1

    mol = Molecule(new_id, name)
    push!(s.molecules, [mol.id, mol.name])

    return mol
end

function add_atom!(s::System, atom::Atom)
    new_id = isempty(s.atoms) ? 1 : maximum(s.atoms.id)+1

    atom = (; atom..., id=new_id)

    push!(s.atoms, atom)
end