import MolecularGraph: 
    GraphMol, graphmol, todict,
    SDFileAtom, SDFileBond, 
    SmilesAtom, SmilesBond

function _atom_to_molgraph(a)
    Dict{String, Any}(
        "stereo" => :unspecified, # TODO: handle stereo chemistry
        "multiplicity" => 1,      # TODO: handle multiplicities
        "symbol" => a.element,
        "charge" => has_property(
            copy(a), "charge") 
                ? get_property(copy(a), "charge") 
                : 0,              # TODO: find a better way to map DataFrameRow to Atom
        "mass"   => nothing,      # TODO: handle masses
        "coords" => a.r
    )
end

function _bond_to_molgraph_edge(b)
    (b.a1, b.a2)
end

function _bond_to_molgraph_attr(b)
    Dict{String, Any}(
        "notation" => 0,            # TODO: handle notation
        "stereo"   => :unspecified, # TODO: handle stereo chemistry
        "order"    => Int(b.order)
    )
end

function Base.convert(
        ::Type{GraphMol{SDFileAtom, SDFileBond}},
        mol::AbstractMolecule{T}
    ) where {T<:Real}

    # create an intermediate dictionary to map from our data structures
    # to those of MolecularGraph
    d = Dict{String, Any}(
        "nodetype"   => "SDFileAtom",
        "nodeattrs"  => map(_atom_to_molgraph, eachrow(mol.atoms)),
        "edgetype"   => "SDFileBond",
        "edges"      => map(_bond_to_molgraph_edge, eachrow(mol.bonds)),
        "edgeattrs"  => map(_bond_to_molgraph_attr, eachrow(mol.bonds)),
        "cache"      => Dict{Any, Any}(),
        "attributes" => mol.properties
    )

    graphmol(d)
end

function _molgraph_to_atom((i, a)::Tuple{Int, SDFileAtom}, T)
    (
        number  = i,
        name    = "$(a.symbol)$(i)",
        element = getproperty(Elements, a.symbol),
        atomtype = "",
        r = Vector3{T}(a.coords),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        has_velocity = false,
        has_force = false,
        frame_id = 1,
        properties = Properties(
            "charge"       => a.charge,
            "multiplicity" => a.multiplicity,
            "mass"         => a.mass,
            "stereo"       => a.stereo
        )
    )
end

function _molgraph_to_atom((i, a)::Tuple{Int, SmilesAtom}, T)
    (
        number  = i,
        name    = "$(a.symbol)$(i)",
        element = getproperty(Elements, a.symbol),
        atomtype = "",
        r = zeros(Vector3{T}),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        has_velocity = false,
        has_force = false,
        frame_id = 1,
        properties = Properties(
            "charge"       => a.charge,
            "multiplicity" => a.multiplicity,
            "mass"         => a.mass,
            "stereo"       => a.stereo,
            "is_aromatic"  => a.isaromatic
        )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SDFileBond}})
    (
        a1 = e[1],
        a2 = e[2],
        order = b.order,
        properties = Properties(
            "notation" => b.notation,
            "stereo"   => b.stereo,
        )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SmilesBond}})
    (
        a1 = e[1],
        a2 = e[2],
        order = b.order,
        properties = Properties(
            "is_aromatic" => b.isaromatic,
            "direction"   => b.direction,
            "stereo"      => b.stereo
        )
    )
end

function Base.convert(
    ::Type{Molecule{T}},
    mol::GraphMol{GMAtom, GMBond}) where {T<:Real, GMAtom, GMBond}
    
    d = todict(mol)
    
    Molecule{T}("",
        DataFrame((t -> _molgraph_to_atom(t, T)).(enumerate(mol.nodeattrs))),
        DataFrame(_molgraph_to_bond.(enumerate(zip(mol.edges, mol.edgeattrs)))),
        Dict(string(k) => v for (k, v) in mol.attributes)
    )
end