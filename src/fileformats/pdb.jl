using BioStructures: 
    read,
    writepdb,
    collectatoms, 
    collectchains, 
    collectresidues, 
    PDB,
    ProteinStructure,
    Model,
    unsafe_addatomtomodel!,
    AtomRecord,
    fixlists!

using Printf

export load_pdb, write_pdb, is_hetero_atom

function is_hetero_atom(a::Atom{T}) where {T<:Real}
    f = parent_fragment(a)

    if !isnothing(f) && is_amino_acid(f)
        return false
    end

    true
end 

function parse_element_string(es::String)
    result = Elements.Unknown

    # handle special cases
    if es == "D"
        result = Elements.H
    elseif es == "X"
        result = Elements.Unknown
    else
        try
            result = parse(Elements, titlecase(es))
        catch e
            @warn "BiochemicalAlgorithms::PDB::parse_element_string: could not parse element from $(es); returning Unknown"
        end
    end

    return result
end

function extract_element(pdb_element::String, atom_name::String)
    element = Elements.Unknown

    # if the element is non-empty, it takes precedence
    if !isempty(strip(pdb_element))
        element = parse_element_string(pdb_element)
    end

    if element == Elements.Unknown && pdb_element != "X"
        # this approach is taken from the original BALL PDB parser

        # try to reconstruct the element from the atom name
        # NOTE: this leads to wrong results if names not compatible with the PDB standard are used, 
        #       such as HE12, which would be interpreted as He
        if (atom_name[1] == ' ' || isdigit(atom_name[1]))
            if atom_name[2] == ' '
                @warn "BiochemicalAlgorithms::PDB::extract_element: could not parse element from atomname $(atom_name); returning Unknown"
                element = Elements.Unknown
            else
                element = parse_element_string(string(atom_name[2]))
            end
        else
            element = parse_element_string(atom_name[1:2])
        end
    end

    return element
end

function Base.convert(::Type{System{T}}, orig_pdb::ProteinStructure) where {T<:Real}
    orig_df  = DataFrame(collectatoms(orig_pdb))

    # then, convert to our representation
    sys = System{T}(orig_pdb.name)
    mol = Molecule(sys, sys.name)

    ### convert the atom positions
    r = Vector3{T}.(T.(orig_df.x), T.(orig_df.y), T.(orig_df.z))

    ### extracting the elements is a little more complicated than it could be,
    ### as the DataFrame-conversion strips whitespaces from atom names
    elements = extract_element.(orig_df.element, getproperty.(collectatoms(orig_pdb, expand_disordered=true), :name))

    atoms = DataFrame(
        number=orig_df.serial, 
        name=orig_df.atomname,
        element=elements,
        r = r,
        formal_charge = orig_df.charge
    )

    # convert other columns of interest to atom properties
    atoms.properties = Properties.(
        collect(
                zip(
                    Pair.(:tempfactor,            orig_df.tempfactor),
                    Pair.(:occupancy,             orig_df.occupancy),
                    Pair.(:is_hetero_atom,        orig_df.ishetero),
                    Pair.(:insertion_code,        orig_df.inscode),
                    Pair.(:alternate_location_id, orig_df.altlocid)
                )
        )
    )

    atoms.frame_id = orig_df.modelnumber
    atoms.chain_id = orig_df.chainid
    atoms.fragment_id = orig_df.resnumber
    atoms.inscode = orig_df.inscode

    # note: we will remove this column as soon as we have filtered out alternates
    atoms.altlocid = orig_df.altlocid

    # collect fragment information
    for orig_chain in collectchains(orig_pdb)
        chain = Chain(mol, orig_chain.id)
        for orig_frag in collectresidues(orig_chain)
            Fragment(chain, orig_frag.number, orig_frag.name, Properties([
                :is_hetero_fragment => orig_frag.het_res,
                :insertion_code => orig_frag.ins_code
            ]))
        end
    end

    # now, handle alternate location ids

    # we try to be as tolerant as possible, as the PDB file format does not
    # seem to formalize many restrictions here.
    # general idea:
    #   - for each atom that has alternative locations, find them
    #   - find the smallest alternate location id and use this as the base case
    #   - store all other variants as properties
    all_altlocs = groupby(filter(:altlocid => !=(' '), atoms, view=true), [:fragment_id, :name])
    for altlocs in all_altlocs
        sorted_altlocs = sort(altlocs, :altlocid, view=true)

        base_case = sorted_altlocs[1, :]
        base_case.altlocid = ' '

        if nrow(sorted_altlocs) > 1
            for altloc in eachrow(sorted_altlocs[2:end, :])
                atoms.properties[base_case.number][Symbol("alternate_location_$(altloc.altlocid)")] = altloc
            end
        end
    end

    # drop all alternates
    atoms = filter(:altlocid => ==(' '), atoms)

    # add all remaining atoms to the system
    grp_atoms = groupby(atoms, [:chain_id, :fragment_id, :inscode])
    for frag in eachfragment(mol)
        for atom in eachrow(grp_atoms[(
            chain_id = parent_chain(frag).name,
            fragment_id = frag.number,
            inscode = frag.properties[:insertion_code]
        )])
            push!(frag, AtomTuple{T}(atom.number, atom.element;
                name = atom.name,
                r = atom.r,
                properties = atom.properties,
            ), frame_id = atom.frame_id)
        end
    end

    sys
end

# Note: models are stored as frames
# TODO: how to handle disordered atoms properly?
function load_pdb(fname::String, T=Float32)
    # first, read the structure using BioStructures.jl
    orig_pdb = read(fname, PDB)
    convert(System{T}, orig_pdb)
end

function _to_atom_record(a::Atom{T}) where {T<:Real}
    # TODO: handle alternative location identifiers!
    f = parent_fragment(a)

    AtomRecord(
        is_hetero_atom(a),
        a.number,
        @sprintf("%4s", a.name),
        ' ',
        f.name,
        parent_chain(a).name,
        f.number,
        get_property(a, :insertion_code, ' '),
        Vector{Float64}(a.r),
        get_property(a, :occupancy, 1.0),
        get_property(a, :tempfactor, 0.0),
        string(a.element),
        string(a.formal_charge)
    )
end

function Base.convert(::Type{ProteinStructure}, ac::AbstractAtomContainer{T}) where {T<:Real}
    # Build a ProteinStructure and add to it incrementally
    struc = ProteinStructure(ac.name)

    # figure out if all molecules in the atom container have the same frames
    sys_frame_ids = frame_ids(ac)

    if ac isa System{T}
        for m in molecules(ac)
            if frame_ids(m) != sys_frame_ids
                error("pdb.jl: cannot convert System containing molecules with different numbers of frames")
            end
        end
    end

    for (i, frame_id) in enumerate(sys_frame_ids)
        struc[i] = Model(i, struc)
       
        for a in atoms(ac; frame_id=frame_id)
            unsafe_addatomtomodel!(
                struc[i],
                _to_atom_record(a)
            )
        end
    end

    fixlists!(struc)
    struc
end

function write_pdb(fname::String, ac::AbstractAtomContainer{T}) where {T<:Real}
    ps = convert(ProteinStructure, ac)
    writepdb(fname, ps)
end