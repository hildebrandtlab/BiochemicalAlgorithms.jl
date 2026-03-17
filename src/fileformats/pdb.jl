export
    is_hetero_atom,
    load_mmcif,
    load_pdb,
    write_mmcif,
    write_pdb

using Printf

function is_hetero_atom(a::Atom)
    f = parent_fragment(a)
    has_flag(a, :is_hetero_atom) || isnothing(f) || !is_amino_acid(f)
end

function parse_element_string(es::AbstractString)
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


"""
    load_pdb(io::IO, ::Type{T} = Float32) -> System{T}
    load_pdb(filename::AbstractString, ::Type{T} = Float32) -> System{T}

Read a PDB file.

!!! note
    Models are stored as frames, using the model number as `frame_id`.
"""
function load_pdb(
        fname_io::Union{AbstractString, IO},
        ::Type{T} = Float32;
        keep_metadata::Bool=true,
        strict_line_checking::Bool=false,
        selected_model::Int=-1,
        ignore_xplor_pseudo_atoms::Bool=true,
        create_coils::Bool=true) where {T <: Real}

    pdblines = readlines(fname_io)

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

    if pdb_info.alternate_location_warning
        @warn "load_pdb: alternate locations other than A are currently not supported. Affected records have been ignored!"
    end

    sys
end

"""
    load_mmcif(io::IO, ::Type{T} = Float32) -> System{T}
    load_mmcif(filename::AbstractString, ::Type{T} = Float32) -> System{T}

Read a PDBx/mmCIF file.

!!! note
    Models are stored as frames, using the model number as `frame_id`.
"""
function load_mmcif(fname_io::Union{AbstractString, IO}, ::Type{T} = Float32;
        create_coils::Bool = true) where {T <: Real}
    MMCIFDetails.read_mmcif(fname_io, T; create_coils=create_coils)
end

function write_pdb(io::IO, ac::Union{Chain{T}, Fragment{T}}) where T
        # get the PDBInfo object
    pdb_info = get_property(ac, :PDBInfo, PDBDetails.PDBInfo{T}())
    pdb_info.writer_stats = PDBDetails.PDBWriterStats()

    PDBDetails.write_coordinate_section(io, pdb_info, ac; coordinate_only=true)
    PDBDetails.write_record(io, pdb_info, PDBDetails.RECORD_TAG_END)
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
    write_pdb(io::IO, ac::AbstractAtomContainer)
    write_pdb(filename::AbstractString, ac::AbstractAtomContainer)

Save an atom container as PDB file.
"""
function write_pdb(fname::AbstractString, ac::AbstractAtomContainer)
    open(io -> write_pdb(io, ac), fname, "w")
end

"""
    write_mmcif(io::IO, ac::AbstractAtomContainer)
    write_mmcif(filename::AbstractString, ac::AbstractAtomContainer)

Save an atom container as PDBx/mmCIF file.
"""
function write_mmcif(io::IO, ac::AbstractAtomContainer)
    MMCIFDetails.write_mmcif_impl(io, ac)
end

function write_mmcif(fname::AbstractString, ac::AbstractAtomContainer)
    open(io -> write_mmcif(io, ac), fname, "w")
end
