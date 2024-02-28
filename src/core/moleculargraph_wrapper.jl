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
        "coords" => a.r,
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

function _molgraph_to_atom((i, a)::Tuple{Int, SDFileAtom}, T)
    (
        number  = i,
        element = getproperty(Elements, a.symbol),
        name    = "$(a.symbol)$(i)",
        atomtype = "",
        r = Vector3{T}(a.coords),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        formal_charge = a.charge,
        charge = zero(T),
        radius = zero(T),
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass,
            :stereo       => a.stereo
        )
    )
end

function _molgraph_to_atom((i, a)::Tuple{Int, SmilesAtom}, T)
    (
        number  = i,
        element = getproperty(Elements, a.symbol),
        name    = "$(a.symbol)$(i)",
        atomtype = "",
        r = zeros(Vector3{T}),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        formal_charge = a.charge,
        charge = zero(T),
        radius = zero(T),
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass,
            :stereo       => a.stereo,
            :is_aromatic  => a.isaromatic
        )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SDFileBond}}, mol)
    at = atoms(mol)
    (
        a1 = only(filter(atom -> atom.number == e[1], at)).idx,
        a2 = only(filter(atom -> atom.number == e[2], at)).idx,
        order = b.order,
        properties = Properties(
            :notation => b.notation,
            :stereo   => b.stereo,
        )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SmilesBond}}, mol)
    at = atoms(mol)
    (
        a1 = only(filter(atom -> atom.number == e[1], at)).idx,
        a2 = only(filter(atom -> atom.number == e[2], at)).idx,
        order = b.order,
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
    system=default_system()) where {T<:Real, GMAtom, GMBond}
    
    d = todict(mg)
    
    mol = Molecule(system, "", Dict{Symbol, Any}(Symbol(k) => v for (k, v) in mg.attributes))

    for a in (t -> _molgraph_to_atom(t, T)).(enumerate(mg.nodeattrs))
        Atom(mol, a...)
    end

    for b in _molgraph_to_bond.(enumerate(zip(mg.edges, mg.edgeattrs)), Ref(mol))
        Bond(parent_system(mol), b.a1, b.a2, BondOrderType(b.order), b.properties)
    end
    
    mol
end