export AbstractMolecule, Molecule, molecules, molecules_df, eachmolecule, nmolecules, parent_molecule

abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

struct Molecule{T} <: AbstractMolecule{T}
    sys::System{T}
    row::DataFrameRow
end

function Molecule(sys::System{T}, name::String = "", properties::Properties = Properties()) where T
    idx = _next_idx(sys)
    push!(sys.molecules, (idx = idx, name = name, properties = properties))
    _molecule_by_idx(sys, idx)
end

function Molecule(name::String = "", properties::Properties = Properties())
    Molecule(default_system(), name, properties)
end

function Base.getproperty(mol::Molecule, name::Symbol)
    in(name, fieldnames(MoleculeTuple)) && return getproperty(getfield(mol, :row), name)
    getfield(mol, name)
end

function Base.setproperty!(mol::Molecule, name::Symbol, val)
    in(name, fieldnames(MoleculeTuple)) && return setproperty!(getfield(mol, :row), name, val)
    setfield!(mol, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", mol::Molecule) = show(io, getfield(mol, :row))
@inline Base.show(io::IO, mol::Molecule) = show(io, getfield(mol, :row))

@inline Base.parent(mol::Molecule) = mol.sys
@inline parent_system(mol::Molecule) = parent(mol)

# TODO this should also add all related entities
#function Base.push!(sys::System{T}, mol::Molecule{T}) where T
#    Molecule(sys, mol.name, mol.properties)
#    sys
#end

@inline function _molecule_by_idx(sys::System{T}, idx::Int) where T
    Molecule{T}(sys, DataFrameRow(sys.molecules, findfirst(sys.molecules.idx .== idx), :))
end

@inline function _molecules(sys::System)
    sys.molecules
end

@inline function molecules(sys::System)
    collect(eachmolecule(sys))
end

@inline function molecules_df(sys::System{T}) where T
    SystemDataFrame{T}(sys, view(_molecules(sys), :, :))
end

@inline function eachmolecule(sys::System{T}) where T
    (Molecule{T}(sys, row) for row in eachrow(_molecules(sys)))
end

function nmolecules(sys::System)
    nrow(_molecules(sys))
end

#=
    Molecule atoms
=#
@inline _atoms(mol::Molecule; kwargs...) = _atoms(mol.sys, molecule_id = mol.idx, kwargs...)
@inline atoms(mol::Molecule; kwargs...) = atoms(mol.sys; molecule_id = mol.idx, kwargs...)
@inline atoms_df(mol::Molecule; kwargs...) = atoms_df(mol.sys; molecule_id = mol.idx, kwargs...)
@inline eachatom(mol::Molecule; kwargs...) = eachatom(mol.sys; molecule_id = mol.idx, kwargs...)
@inline natoms(mol::Molecule; kwargs...) = natoms(mol.sys; molecule_id = mol.idx, kwargs...)

@inline function Base.push!(mol::Molecule{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(mol.sys, atom; molecule_id = mol.idx, kwargs...)
    mol
end

#=
    Molecule bonds
=#
@inline _bonds(mol::Molecule; kwargs...) = _bonds(mol.sys; molecule_id = mol.idx, kwargs...)
@inline bonds(mol::Molecule; kwargs...) = bonds(mol.sys; molecule_id = mol.idx, kwargs...)
@inline bonds_df(mol::Molecule; kwargs...) = bonds_df(mol.sys; molecule_id = mol.idx, kwargs...)
@inline eachbond(mol::Molecule; kwargs...) = eachbond(mol.sys; molecule_id = mol.idx, kwargs...)
@inline nbonds(mol::Molecule; kwargs...) = nbonds(mol.sys; molecule_id = mol.idx, kwargs...)

@inline function Base.push!(mol::Molecule, bond::BondTuple)
    push!(mol.sys, bond)
    mol
end
