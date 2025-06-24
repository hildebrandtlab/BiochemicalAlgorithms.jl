export
    SMARTSQuery

struct SMARTSQuery
    query::String
    query_graph::MolecularGraph.SMARTSMolGraph

    function SMARTSQuery(query::AbstractString)
        new(query, MolecularGraph.smartstomol(query))
    end
end

@inline function _to_substructure(name::AbstractString, mol::AbstractAtomContainer, idxset::Set{Int}; adjacent_bonds::Bool = false)
    filter_atoms(atom -> atom.idx ∈ idxset, mol; name, adjacent_bonds)
end

function Base.match(query::SMARTSQuery, mol::AbstractAtomContainer; adjacent_bonds::Bool = false)
    mg_mol = convert(MolecularGraph.SDFMolGraph, mol)

    a2idx  = Dict(mg_mol.vprops[a] => idx for (a, idx) in _molgraph_atom_to_idx(mg_mol))
    matches = MolecularGraph.substruct_matches(mg_mol, query.query_graph)

    [_to_substructure("$(query.query) on $(mol.name)", mol, Set(a2idx[mg_mol.vprops[a]] for (a, _) in m); adjacent_bonds) for m in matches]
end
