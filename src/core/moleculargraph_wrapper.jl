import MolecularGraph: 
    GraphMol, graphmol, todict,
    SDFileAtom, SDFileBond, 
    SmilesAtom, SmilesBond

function _atom_to_molgraph(a)
    Dict{String, Any}(
        "stereo" => :unspecified, # TODO: handle stereo chemistry
        "multiplicity" => 1,      # TODO: handle multiplicities
        "symbol" => a.element,
        "charge" => a.formal_charge,
        "mass"   => nothing,      # TODO: handle masses
        "coords" => ustrip.(a.r),
        "idx"    => a.idx
    )
end

function _bond_to_molgraph_edge(b, idx_to_molgraph_atom)
    (idx_to_molgraph_atom[b.a1], idx_to_molgraph_atom[b.a2])
end

function _bond_to_molgraph_attr(b)
    Dict{String, Any}(
        "notation" => 0,            # TODO: handle notation
        "stereo"   => :unspecified, # TODO: handle stereo chemistry
        "order"    => b.order <= BondOrder.Quadruple ? Int(b.order) : 1, # TODO: handle aromatic bonds correctly
        "idx"      => b.idx
    )
end

function Base.convert(
        ::Type{GraphMol{SDFileAtom, SDFileBond}},
        mol::AbstractAtomContainer{T}
    ) where {T<:Real}

    # create an intermediate dictionary to map from our data structures
    # to those of MolecularGraph
    molgraph_atoms = map(_atom_to_molgraph, atoms(mol))
    idx_to_molgraph_atom = Dict(
        a["idx"] => i for (i,a) in enumerate(molgraph_atoms)
    )
    molgraph_atom_to_idx = Dict(
        v => k for (k,v) in idx_to_molgraph_atom
    )

    d = Dict{String, Any}(
        "nodetype"   => "SDFileAtom",
        "nodeattrs"  => molgraph_atoms,
        "edgetype"   => "SDFileBond",
        "edges"      => map(b -> _bond_to_molgraph_edge(b, idx_to_molgraph_atom), bonds(mol)),
        "edgeattrs"  => map(_bond_to_molgraph_attr, bonds(mol)),
        "cache"      => Dict{Any, Any}(),
        "attributes" => mol.properties ∪ Dict("atom_idx" => molgraph_atom_to_idx)
    )

    graphmol(d)
end

@inline function _molgraph_to_atom(mol::Molecule{T}, (i, a)::Tuple{Int, SDFileAtom}) where T
    Atom(mol, i, getproperty(Elements, a.symbol);
        name = "$(a.symbol)$(i)",
        r = Vector3{T}(a.coords) * u"Å",
        formal_charge = a.charge,
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass,
            :stereo       => a.stereo
        )
    )
end

@inline function _molgraph_to_atom(mol::Molecule{T}, (i, a)::Tuple{Int, SmilesAtom}) where T
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

@inline function _molgraph_to_bond(mol::Molecule, (i, (e, b))::Tuple{Int, Tuple{Any, SDFileBond}})
    at = atoms(mol)
    Bond(
        mol,
        only(filter(atom -> atom.number == e[1], at)).idx,
        only(filter(atom -> atom.number == e[2], at)).idx,
        BondOrderType(b.order);
        properties = Properties(
            :notation => b.notation,
            :stereo   => b.stereo,
        )
    )
end

@inline function _molgraph_to_bond(mol::Molecule, (i, (e, b))::Tuple{Int, Tuple{Any, SmilesBond}})
    at = atoms(mol)
    Bond(
        mol,
        only(filter(atom -> atom.number == e[1], at)).idx,
        only(filter(atom -> atom.number == e[2], at)).idx,
        BondOrderType(b.order);
        properties = Properties(
            :is_aromatic => b.isaromatic,
            :direction   => b.direction,
            :stereo      => b.stereo
        )
    )
end

function Base.convert(
    ::Type{Molecule{T}},
    mg::GraphMol{GMAtom, GMBond};
    system=default_system()
) where {T, GMAtom, GMBond}
    mol = Molecule(system; properties = Properties(Symbol(k) => v for (k, v) in mg.attributes))
    _molgraph_to_atom.(Ref(mol), enumerate(mg.nodeattrs))
    _molgraph_to_bond.(Ref(mol), enumerate(zip(mg.edges, mg.edgeattrs)))
    mol
end
