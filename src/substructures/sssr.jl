import MolecularGraph: 
    GraphMol, 
    SDFileAtom, 
    SDFileBond,
    sssr,
    isringatom

export find_sssr, is_ring_atom

function _filter_bonds(ac::AbstractAtomContainer{T}) where {T<:Real}
     # filter out cystein bridges and h bridges (TODO!)
     new_atoms = Tables.materializer(AtomTable{T})(_atoms(ac))
     new_bonds = Tables.materializer(BondTable)(
        TableOperations.filter(row ->
            !get(row.properties, "DISULPHIDE_BOND", false),
            _bonds(ac)
        )
     )

    convert(GraphMol{SDFileAtom, SDFileBond}, Substructure(ac.name, ac, new_atoms, new_bonds, ac.properties))
end

function find_sssr(ac::AbstractAtomContainer{T}) where {T<:Real}
    mg = _filter_bonds(ac)

    mg_sssr = sssr(mg)

    map(r->map(a->atom_by_idx(ac isa System ? ac : parent_system(ac), mg.attributes[:atom_idx][a]), r), mg_sssr)
end

function is_ring_atom(ac::AbstractAtomContainer{T}) where {T<:Real}
    mg = _filter_bonds(ac)

    isringatom(mg)
end