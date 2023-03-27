export Substructure, filter_atoms

@auto_hash_equals struct Substructure{T<:Real} <: AbstractMolecule{T}
    name::String

    parent::AbstractAtomContainer{T}
    atoms::SubDataFrame
    bonds::SubDataFrame
    properties::Properties

    function Substructure{T}(name,
                             parent, 
                             atoms, 
                             bonds, 
                             properties = Properties()) where {T<:Real}
        new(name, parent, atoms, bonds, properties)
    end
end

Substructure(name,
             parent,
             atoms,
             bonds,
             properties = Properties()) = 
                Substructure{Float32}(name, parent, atoms, bonds, properties)

function filter_atoms(fn, mol; name="", adjacent_bonds=false)
    atom_view = filter(fn, getfield(atoms_df(mol), :df), view=true)
    bond_view = filter(
        [:a1, :a2] => (a1, a2) -> 
            adjacent_bonds ? a1 ∈ atom_view.idx || a2 ∈ atom_view.idx
                           : a1 ∈ atom_view.idx && a2 ∈ atom_view.idx,
        getfield(bonds_df(mol), :df), view=true)

    Substructure(name, mol, atom_view, bond_view)
end

function atoms(ac::Substructure{T}) where {T<:Real}
    collect(eachrow(ac.atoms))
end

function bonds(ac::Substructure{T}) where {T<:Real}
    collect(eachrow(ac.bonds))
end