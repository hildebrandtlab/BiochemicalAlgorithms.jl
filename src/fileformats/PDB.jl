using BioStructures: read, collectatoms, PDB

export load_pdb

function parse_element_string(es::String)
    result = Elements.Unknown

    # handle special cases
    if es == "D"
        result = Elements.H
    elseif es == "X"
        result = Elements.Unknown
    else
        try
            result = parse(Elements, es)
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

# Note: models are stored as frames
# TODO: how to handle disordered atoms properly?
function load_pdb(fname::String, T=Float32)
    # first, read the structure using BioStructures.jl
    orig_pdb = read(fname, PDB)
    orig_df  = DataFrame(collectatoms(orig_pdb))

    # then, convert to our representation

    ### convert the atom positions
    r = Vector3{T}.(T.(orig_df.x), T.(orig_df.y), T.(orig_df.z))

    ### extracting the elements is a little more complicated than it could be,
    ### as the DataFrame-conversion strips whitespaces from atom names
    elements = extract_element.(orig_df.element, getproperty.(collectatoms(orig_pdb, expand_disordered=true), :name))

    atoms = DataFrame(
        number=orig_df.serial, 
        name=orig_df.atomname,
        element=elements
    )

    insertcols!(atoms, :atomtype => "")
    atoms.r = r
    atoms.v .= Ref(Vector3{T}(0.0, 0.0, 0.0))
    atoms.F .= Ref(Vector3{T}(0.0, 0.0, 0.0))
    atoms.has_velocity .= Ref(false)
    atoms.has_force .= Ref(false)
    atoms.frame_id = orig_df.modelnumber

    atoms.residue_id = orig_df.resnumber
    atoms.residue_name= orig_df.resname
    atoms.chain = orig_df.chainid

    # note: we will remove this column as soon as we have filtered out alternates
    atoms.altlocid = orig_df.altlocid

    chains = DataFrame(Chain.(Ref.(unique(orig_df.chainid))))

    # convert other columns of interest to atom properties
    atom_properties = Dict(
        Pair.(orig_df.serial, 
            AtomProperties.(
                collect(
                    zip(
                        Pair.("tempfactor",            orig_df.tempfactor),
                        Pair.("occupancy",             orig_df.occupancy),
                        Pair.("is_hetero_atom",        orig_df.ishetero),
                        Pair.("insertion_code",        orig_df.inscode),
                        Pair.("alternate_location_id", orig_df.altlocid)
                    )
                )
            )
        )
    )

    bonds = DataFrame(Bond[])
    residues = DataFrame(Residue[])

    # now, handle alternate location ids

    # we try to be as tolerant as possible, as the PDB file format does not
    # seem to formalize many restrictions here.
    # general idea:
    #   - for each atom that has alternative locations, find them
    #   - find the smallest alternate location id and use this as the base case
    #   - store all other variants as properties
    all_altlocs = groupby(filter(:altlocid => !=(' '), atoms, view=true), [:residue_id, :name])
    for altlocs in all_altlocs
        sorted_altlocs = sort(altlocs, :altlocid, view=true)

        base_case = sorted_altlocs[1, :]
        base_case.altlocid = ' '

        if nrow(sorted_altlocs) > 1
            for altloc in eachrow(sorted_altlocs[2:end, :])
                atom_properties[base_case.number]["alternate_location_$(altloc.altlocid)"] = altloc
            end
        end
    end

    # drop all alternates
    atoms = filter(:altlocid => ==(' '), atoms)

    # drop the altlocid-column
    select!(atoms, Not(:altlocid))
    
    p = Protein{T}(orig_pdb.name, atoms, bonds, atom_properties, residues, chains)
end
