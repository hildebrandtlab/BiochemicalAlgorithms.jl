export
    find_sssr,
    is_ring_atom

function _filter_bonds(ac::AbstractAtomContainer)
     # filter out cystein bridges and h bridges (TODO!)
     new_atoms = atoms(ac)
     new_bonds = filter(row -> !has_flag(row, :TYPE__DISULPHIDE_BOND), bonds(ac))

    convert(MolecularGraph.SDFMolGraph, Substructure(ac.name, ac, new_atoms, new_bonds, ac.properties))
end

function find_sssr(ac::AbstractAtomContainer)
    mg = _filter_bonds(ac)

    mg_sssr = MolecularGraph.sssr(mg)

    map(r->map(a->atom_by_idx(ac isa System ? ac : parent_system(ac), mg.gprops[:atom_idx][a]), r), mg_sssr)
end

@inline function is_ring_atom(ac::AbstractAtomContainer)
    MolecularGraph.is_in_ring(_filter_bonds(ac))
end
