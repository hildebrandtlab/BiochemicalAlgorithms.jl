@inline function _atom_to_molgraph(a::Atom)
    MolecularGraph.SDFAtom(Dict{String, Any}(
        "multiplicity" => 1,      # TODO: handle multiplicities
        "symbol" => a.element,
        "charge" => a.formal_charge,
        "mass"   => nothing,      # TODO: handle masses
        "coords" => a.r
    ))
end

@inline function _bond_to_molgraph_edge(b::Bond, idx_to_molgraph_atom::Dict{Int,Int})
    # GraphMol assumes src < dst
    e1 = idx_to_molgraph_atom[b.a1]
    e2 = idx_to_molgraph_atom[b.a2]
    e1 < e2 && return Edge(e1, e2)
    Edge(e2, e1)
end

@inline function _bond_to_molgraph_attr(b::Bond)
    MolecularGraph.SDFBond(Dict{String, Any}(
        "isordered" => has_property(b, :is_odered) ? Bool(get_property(b, :is_ordered)) : true,
        "notation" => 0,                                                # TODO: handle notation
        "order"    => b.order <= BondOrder.Quadruple ? Int(b.order) : 1 # TODO: handle aromatic bonds correctly
    ))
end

function Base.convert(
        ::Type{MolecularGraph.SDFMolGraph},
        mol::AbstractAtomContainer{T}
    ) where T

    # create an intermediate dictionary to map from our data structures
    # to those of MolecularGraph
    at = atoms(mol)
    molgraph_atoms = sort!(OrderedDict(a.number => _atom_to_molgraph(a) for a in at))
    idx_to_molgraph_atom = Dict(a.idx => a.number for a in at)
    molgraph_atom_to_idx = Dict(v => k for (k,v) in idx_to_molgraph_atom)

    edges  = collect(Edge{Int}, map(b -> _bond_to_molgraph_edge(b, idx_to_molgraph_atom), bonds(mol)))
    vprops = collect(MolecularGraph.SDFAtom,   values(molgraph_atoms))
    eprops = collect(MolecularGraph.SDFBond,   map(_bond_to_molgraph_attr, bonds(mol)))
    
    MolecularGraph.SDFMolGraph(edges, vprops, eprops;
       gprop_map = merge(mol.properties, Dict(:atom_idx => molgraph_atom_to_idx))
    )
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
    mol = Molecule(system; properties = Properties(mg.gprops))
    foreach(t -> _molgraph_to_atom(mol, t), sort!(OrderedDict(mg.vprops)))
    foreach(t -> _molgraph_to_bond(mol, t), sort!(OrderedDict(mg.eprops)))
    mol
end
