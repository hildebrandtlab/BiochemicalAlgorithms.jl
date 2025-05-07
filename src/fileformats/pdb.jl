export
    is_hetero_atom,
    load_mmcif,
    load_pdb,
    write_mmcif,
    write_pdb

using BioStructures:
    read,
    writemmcif,
    writepdb,
    collectatoms,
    collectchains,
    collectresidues,
    PDBFormat,
    MMCIFFormat,
    MolecularStructure,
    Model,
    unsafe_addatomtomodel!,
    AtomRecord,
    fixlists!

using Printf

function is_hetero_atom(a::Atom)
    f = parent_fragment(a)
    has_flag(a, :is_hetero_atom) || isnothing(f) || !is_amino_acid(f)
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
            result = parse(ElementType, es)
        catch
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

function Base.convert(::Type{System{T}}, orig_pdb::MolecularStructure) where T
    orig_df  = DataFrame(collectatoms(orig_pdb))

    # then, convert to our representation
    sys = System{T}(orig_pdb.name)
    mol = Molecule(sys; name = sys.name)

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

    atoms.flags = map(
        h -> h ? Flags([:is_hetero_atom]) : Flags(),
        orig_df.ishetero
    ) .∪ map(
        d -> d ? Flags([:is_deuterium]) : Flags(),
        orig_df.element .== "D"
    )

    atoms.frame_id = orig_df.modelnumber
    atoms.chain_idx = orig_df.chainid
    atoms.fragment_idx = orig_df.resnumber
    atoms.inscode = orig_df.inscode

    # note: we will remove this column as soon as we have filtered out alternates
    atoms.altlocid = orig_df.altlocid

    # collect fragment information
    for orig_chain in collectchains(orig_pdb)
        chain = Chain(mol; name = orig_chain.id)
        for orig_frag in collectresidues(orig_chain)
            Fragment(chain, orig_frag.number;
                name = orig_frag.name,
                properties = Properties([
                    :is_hetero_fragment => orig_frag.het_res,
                    :insertion_code => orig_frag.ins_code
                ]),
                variant = _fragment_variant(strip(orig_frag.name))
            )
        end
    end

    # now, handle alternate location ids

    # we try to be as tolerant as possible, as the PDB file format does not
    # seem to formalize many restrictions here.
    # general idea:
    #   - for each atom that has alternative locations, find them
    #   - find the smallest alternate location id and use this as the base case
    #   - store all other variants as properties
    all_altlocs = groupby(filter(:altlocid => !=(' '), atoms, view=true), [:chain_idx, :fragment_idx, :name])
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
    grp_atoms = groupby(atoms, [:chain_idx, :fragment_idx, :inscode])
    for frag in fragments(mol)
        for atom in eachrow(grp_atoms[(
            chain_idx = parent_chain(frag).name,
            fragment_idx = frag.number,
            inscode = frag.properties[:insertion_code]
        )])
            Atom(frag, atom.number, atom.element;
                name = atom.name,
                r = atom.r,
                properties = atom.properties,
                flags = atom.flags,
                frame_id = atom.frame_id
            )
        end
    end

    sys
end


"""
    load_pdb(fname::AbstractString, ::Type{T} = Float32) -> System{T}

Read a PDB file.

!!! note
    Models are stored as frames, using the model number as `frame_id`.
"""
function load_pdb(
        filename::AbstractString, 
        ::Type{T} = Float32;
        keep_metadata::Bool=true,
        strict_line_checking::Bool=false,
        selected_model::Int=-1, 
        ignore_xplor_pseudo_atoms::Bool=true,
        create_coils::Bool=true) where {T <: Real}

    pdblines = readlines(filename)

    sys = System{T}()
    mol = Molecule(sys)

    pdb_info = PDBDetails.PDBInfo{T}(selected_model)

    for pl in pdblines
        PDBDetails.handle_record(pl, sys, pdb_info; 
            strict_line_checking=strict_line_checking,
            ignore_xplor_pseudo_atoms=ignore_xplor_pseudo_atoms)
    end

    sys.name = strip(pdb_info.name)
    mol.name = sys.name

    fragment_cache = Dict{PDBDetails.UniqueResidueID, Fragment{T}}()

    for f in fragments(sys)
        fragment_cache[PDBDetails.UniqueResidueID(f.name, parent_chain(f).name, f.number, get_property(f, :insertion_code, " "))] = f
    end

    PDBDetails.postprocess_ssbonds_!(sys, pdb_info, fragment_cache)
    PDBDetails.postprocess_secondary_structures_!(sys, pdb_info, fragment_cache, create_coils)

    if keep_metadata
        # clean up a little... (if we don't do this, comparing systems will lead to a stack overflow)
        pdb_info.current_chain = nothing
        pdb_info.current_residue = nothing

        set_property!(sys, :PDBInfo, pdb_info)
    end

    sys
end

"""
    load_mmcif(fname::AbstractString, ::Type{T} = Float32) -> System{T}

Read a PDBx/mmCIF file.

!!! note
    Models are stored as frames, using the model number as `frame_id`.
"""
function load_mmcif(fname::AbstractString, ::Type{T} = Float32) where {T <: Real}
    # TODO: how to handle disordered atoms properly?

    # first, read the structure using BioStructures.jl
    orig_mmcif = read(fname, MMCIFFormat)
    convert(System{T}, orig_mmcif)
end

function _to_atom_record(a::Atom)
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

function Base.convert(::Type{MolecularStructure}, ac::AbstractAtomContainer{T}) where T
    # Build a MolecularStructure and add to it incrementally
    struc = MolecularStructure(ac.name)

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


function write_pdb(io::IO, ac::AbstractAtomContainer{T}) where T
    # get the PDBInfo object
    pdb_info = get_property(ac, :PDBInfo, PDBDetails.PDBInfo{T}())

    pdb_info.writer_stats = PDBDetails.PDBWriterStats()

    PDBDetails.write_title_section(io, pdb_info)
    PDBDetails.write_primary_structure_section(io, pdb_info, ac)
    PDBDetails.write_heterogen_section(io, pdb_info)
    PDBDetails.write_secondary_structure_section(io, pdb_info, ac)
    PDBDetails.write_connectivity_annotation_section(io, pdb_info, ac)
    PDBDetails.write_miscellaneous_features_section(io, pdb_info)
    PDBDetails.write_crystallographic_section(io, pdb_info)
    PDBDetails.write_coordinate_section(io, pdb_info, ac)
    PDBDetails.write_connectivity_section(io, pdb_info, ac)
    PDBDetails.write_bookkeeping_section(io, pdb_info, ac)
    PDBDetails.write_record(io, pdb_info, PDBDetails.RECORD_TAG_END)
end

"""
    write_pdb(fname::AbstractString, ac::AbstractAtomContainer)

Save an atom container as PDB file.
"""
function write_pdb(fname::AbstractString, ac::AbstractAtomContainer)
    open(io -> write_pdb(io, ac), fname, "w")
end

"""
    write_mmcif(fname::AbstractString, ac::AbstractAtomContainer)

Save an atom container as PDBx/mmCIF file.
"""
function write_mmcif(fname::AbstractString, ac::AbstractAtomContainer)
    ps = convert(MolecularStructure, ac)
    writemmcif(fname, ps)
end
