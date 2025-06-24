@inline function _atom_to_molgraph(a::Atom)
    MolecularGraph.SDFAtom(
        symbol = Symbol(a.element),
        charge = a.formal_charge,
        multiplicity = 1,    # TODO: handle multiplicities
        mass = nothing,      # TODO: handle masses
        coords = collect(a.r)
    )
end

@inline function _bond_to_molgraph_edge(b::Bond, idx_to_molgraph_atom::Dict{Int,Int})
    # MolecularGraph assumes src < dst (cf. MolecularGraph.u_edge)
    e1 = idx_to_molgraph_atom[b.a1]
    e2 = idx_to_molgraph_atom[b.a2]
    e1 < e2 ? Edge(e1, e2) : Edge(e2, e1)
end

@inline function _bond_to_molgraph_attr(b::Bond)
    MolecularGraph.SDFBond(
        isordered = has_property(b, :is_odered) ? Bool(get_property(b, :is_ordered)) : true,
        notation = 0,                                             # TODO: handle notation
        order = b.order <= BondOrder.Quadruple ? Int(b.order) : 1 # TODO: handle aromatic bonds correctly
    )
end

function Base.convert(
    ::Type{MolecularGraph.SDFMolGraph},
    mol::AbstractAtomContainer
)
    at = atoms(mol)

    # create an intermediate dictionary to map from our data structures to those of MolecularGraph
    # NOTE: the number field of the atoms is not guaranteed to be unique and the atom numbering
    # is not guaranteed to be 1-based and contiguous. So, we re-number the atoms here in order.
    molgraph_atoms = OrderedDict(ai => _atom_to_molgraph(a) for (ai, a) in enumerate(at))
    idx_to_molgraph_atom = Dict(a.idx => ai for (ai, a) in enumerate(at))
    molgraph_atom_to_idx = Dict(v => k for (k,v) in idx_to_molgraph_atom)

    edges    = collect(Edge{Int}, map(b -> _bond_to_molgraph_edge(b, idx_to_molgraph_atom), bonds(mol)))
    sdfatoms = collect(MolecularGraph.SDFAtom,   values(molgraph_atoms))
    sdfbonds = collect(MolecularGraph.SDFBond,   map(_bond_to_molgraph_attr, bonds(mol)))

    mg = MolecularGraph.SDFMolGraph(edges, sdfatoms, sdfbonds;
        on_init = MolecularGraph.sdf_on_init!,
        on_update = MolecularGraph.sdf_on_update!
    )
    mg["atom_idx"] = JSON3.write(molgraph_atom_to_idx)
    mg
end

@inline function _molgraph_atom_to_idx(mg::MolecularGraph.MolGraph)
    JSON3.read(mg["atom_idx"], Dict{Int, Int})
end

@inline function _molgraph_to_atom(mol::Molecule{T}, (i, a)::Pair{Int, MolecularGraph.SDFAtom}) where T
    Atom(mol, i, getproperty(Elements, a.symbol);
        name = "$(a.symbol)$(i)",
        r = Vector3{T}(a.coords),
        formal_charge = a.charge,
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass
        )
    )
end

@inline function _molgraph_to_atom(mol::Molecule, (i, a)::Pair{Int, MolecularGraph.SMILESAtom})
    Atom(mol, i, getproperty(Elements, a.symbol);
        name = "$(a.symbol)$(i)",
        formal_charge = a.charge,
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass,
            :stereo       => a.stereo,
            :is_aromatic  => a.isaromatic
        )
    )
end

@inline function _molgraph_to_bond(mol::Molecule, (e, b)::Pair{E, MolecularGraph.SDFBond}) where E
    at = atoms(mol)
    Bond(
        mol,
        only(filter(atom -> atom.number == src(e), at)).idx,
        only(filter(atom -> atom.number == dst(e), at)).idx,
        BondOrderType(b.order);
        properties = Properties(
            :notation => b.notation
        )
    )
end

@inline function _molgraph_to_bond(mol::Molecule, (e, b)::Pair{E, MolecularGraph.SMILESBond}) where E
    at = atoms(mol)
    Bond(
        mol,
        only(filter(atom -> atom.number == src(e), at)).idx,
        only(filter(atom -> atom.number == dst(e), at)).idx,
        BondOrderType(b.order);
        properties = Properties(
            :is_aromatic => b.isaromatic,
            :direction   => b.direction
        )
    )
end

function Base.convert(
    ::Type{Molecule{T}},
    mg::MolecularGraph.MolGraph;
    system = System{T}()
) where T
    mol = Molecule(system; properties = Properties(:metadata => mg.gprops.metadata))
    foreach(t -> _molgraph_to_atom(mol, t), sort!(OrderedDict(mg.vprops)))
    foreach(t -> _molgraph_to_bond(mol, t), sort!(OrderedDict(mg.eprops)))
    mol
end
