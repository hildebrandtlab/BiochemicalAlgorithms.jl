using BioStructures: read, collectatoms, PDB

export load_pdb

function extract_element(pdb_element::String, atom_name::String)
    element = Elements.Unknown

    try
        # if the element is non-empty, it takes precedence
        if !isempty(strip(pdb_element))
            element = parse(Elements, pdb_element)
        end

        # this approach is taken from the original BALL PDB parser

        # try to reconstruct the element from the atom name
        # NOTE: this leads to wrong results if names not compatible with the PDB standard are used, 
        #       such as HE12, which would be interpreted as He
        if (atom_name[1] == ' ' || isdigit(atom_name[0]))
            if atom_name[2] == ' '
                @warn "BiochemicalAlgorithms::PDB::extract_element: could not parse element from atomname $(atom_name); returning Unknown"
                element = Elements.Unknown
            else
                element = parse(Elements, string(atom_name[2]))
            end
        else
            element = parse(Elements, atom_name[1:2])
        end
    catch e
        @warn "BiochemicalAlgorithms::PDB::extract_element: could not parse element from atom $(atom_name), element $(pdb_element); returning Unknown"

        element = Elements.Unknown
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

    bonds = DataFrame(Bond[])
    residues = DataFrame(Residue[])

    p = Protein{T}(orig_pdb.name, atoms, bonds, residues)
end
