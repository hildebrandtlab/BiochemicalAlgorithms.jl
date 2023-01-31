import MolecularGraph: GraphMol, SDFileAtom, SDFileBond, graphmol

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
        mol::AbstractMolecule
    )

    # create an intermediate dictionary to map from our data structures
    # to those of MolecularGraph
    d = Dict{String, Any}(
        "nodetype"   => "SDFileAtom",
        "nodeattrs"  => map(_atom_to_molgraph, eachrow(mol.atoms)),
        "edgetype"   => "SDFileBond",
        "edges"      => map(_bond_to_molgraph_edge, eachrow(mol.bonds)),
        "edgeattrs"  => map(_bond_to_molgraph_attr, eachrow(mol.bonds)),
        "cache"      => Dict{Any, Any}(),
        "attributes" => Dict{Symbol, Any}()
        
    )

    graphmol(d)
end