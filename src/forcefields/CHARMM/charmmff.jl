export
    CharmmFF

get_charmm_default_options(::Type{T} = Float32) where {T <: Real} = Dict{Symbol, Any}(
    :nonbonded_cutoff               => T(20.0),
    :vdw_cutoff                     => T(9.0),
    :vdw_cuton                      => T(7.0),
    :electrostatic_cutoff           => T(9.0),
    :electrostatic_cuton            => T(7.0),
    :solvation_cutoff               => T(9.0),
    :solvation_cuton                => T(7.0),
    :scaling_vdw_1_4                => T(1.0),
    :scaling_electrostatic_1_4      => T(1.0),
    :distance_dependent_dielectric  => true,
    :use_eef1                       => true,
    :remove_extended_hydrogens      => true,
    :assign_charges                 => true,
    :assign_typenames               => true,
    :assign_types                   => true,
    :overwrite_nonzero_charges      => true,
    :overwrite_typenames            => false,
    :periodic_boundary_conditions   => false,
    :periodic_box_width             => T(100.0),
    :periodic_box_height            => T(100.0),
    :periodic_box_depth             => T(100.0),
    :max_number_of_unassigned_atoms => typemax(Int32)
)

# reads the CHARMM nonbonded options (CUTNB, CTOFNB, CTONNB, E14FAC, ATOM/CDIEL/RDIEL)
# stored in the LennardJones section of the parameter file
function _charmm_file_options(params::CharmmFFParameters, ::Type{T}) where T
    file_options = Dict{Symbol, Any}()

    haskey(params.sections, "LennardJones") || return file_options

    props = Dict(uppercase(k) => strip(v) for (k, v) in extract_section(params, "LennardJones").properties)

    if haskey(props, "CUTNB")
        file_options[:nonbonded_cutoff] = parse(T, props["CUTNB"])
    end

    if haskey(props, "CTOFNB")
        cutoff = parse(T, props["CTOFNB"])
        file_options[:vdw_cutoff] = file_options[:electrostatic_cutoff] = file_options[:solvation_cutoff] = cutoff
    end

    if haskey(props, "CTONNB")
        cuton = parse(T, props["CTONNB"])
        file_options[:vdw_cuton] = file_options[:electrostatic_cuton] = file_options[:solvation_cuton] = cuton
    end

    if haskey(props, "E14FAC")
        file_options[:scaling_electrostatic_1_4] = parse(T, props["E14FAC"])
    end

    if get(props, "ATOM", "") == "CDIEL" || get(props, "CDIEL", "") == "true"
        file_options[:distance_dependent_dielectric] = false
    elseif get(props, "ATOM", "") == "RDIEL" || get(props, "RDIEL", "") == "true"
        file_options[:distance_dependent_dielectric] = true
    end

    file_options
end

# CHARMM's united atom parameter sets (e.g., param19) model hydrogens bound to carbon implicitly
# through "extended" atom types; this removes all hydrogens bound to atoms of such types
function _remove_extended_hydrogens!(ff::ForceField{T}) where T
    atom_types = extract_section(ff.parameters, "AtomTypes").data

    extended_types = Set{String}(
        String(r.type) for r in eachrow(atom_types) if r.extended == true || r.extended == "true"
    )

    isempty(extended_types) && return

    templates = ff.atom_type_templates

    function _template_type(atom)
        for name in (
            get_full_name(atom),
            get_full_name(atom, FullNameType.NO_VARIANT_EXTENSIONS),
            "*:" * atom.name
        )
            haskey(templates, name) && return templates[name].type_name
        end
        nothing
    end

    hydrogens = Set{Int}()
    for atom in atoms(ff.system)
        atom.element == Elements.H && continue

        if _template_type(atom) ∈ extended_types
            for bond in bonds(atom)
                partner = get_partner(bond, atom)
                if partner.element == Elements.H
                    push!(hydrogens, partner.idx)
                end
            end
        end
    end

    isempty(hydrogens) && return

    for idx in hydrogens
        delete!(atom_by_idx(ff.system, idx))
    end

    @info "CharmmFF: deleted $(length(hydrogens)) hydrogen atoms (using CHARMM united atoms instead)."
end

"""
    CharmmFF(
        ::AbstractAtomContainer{T},
        param_file::AbstractString = ball_data_path("forcefields/CHARMM/param19_eef1.ini")
    )

Initializes a CHARMM force field for the given atom container and the given parameter
file (default: united-atom PARAM19 with EEF1 implicit solvation). The all-atom parameter
set PARAM22 is available as `ball_data_path("forcefields/CHARMM/param22.ini")`; it uses
CHARMM hydrogen names (e.g., `HN` for the backbone amide hydrogen) and contains no residue
improper torsion definitions, so improper torsions are disabled for PARAM22.

Nonbonded options stored in the parameter file (`CUTNB`, `CTOFNB`, `CTONNB`, `E14FAC`,
`ATOM`) override the defaults below; keyword arguments override both.

!!! note
    For united-atom parameter sets, all hydrogen atoms bound to atoms with an
    `extended` atom type are **deleted from the system** during setup, unless
    `remove_extended_hydrogens = false` is given.

!!! note
    In contrast to [`AmberFF`](@ref), the 1-4 scaling factors are multiplied with
    (not divided by) the respective 1-4 interaction energies.

# Supported keyword arguments
 - `nonbonded_cutoff::T = 20`
 - `vdw_cutoff::T = 9`
 - `vdw_cuton::T = 7`
 - `electrostatic_cutoff::T = 9`
 - `electrostatic_cuton::T = 7`
 - `solvation_cutoff::T = 9`
 - `solvation_cuton::T = 7`
 - `scaling_vdw_1_4::T = 1`
 - `scaling_electrostatic_1_4::T = 1`
 - `distance_dependent_dielectric::Bool = true`
 - `use_eef1::Bool = true`
 - `remove_extended_hydrogens::Bool = true`
 - `assign_charges::Bool = true`
 - `assign_typenames::Bool = true`
 - `assign_types::Bool = true`
 - `overwrite_nonzero_charges::Bool = true`
 - `overwrite_typenames::Bool = false`
 - `periodic_boundary_conditions::Bool = false`
 - `periodic_box_width::T = 100`
 - `periodic_box_height::T = 100`
 - `periodic_box_depth::T = 100`
 - `max_number_of_unassigned_atoms::Int = typemax(Int32)`
"""
function CharmmFF(
    ac::AbstractAtomContainer{T},
    param_file::AbstractString=ball_data_path("forcefields/CHARMM/param19_eef1.ini");
    constrained_atoms=Vector{Int}(),
    kwargs...
) where T

    charmm_params = CharmmFFParameters(param_file, T)
    options = get_charmm_default_options(T)

    merge!(options, _charmm_file_options(charmm_params, T))

    for (key, value) in kwargs
        if !haskey(options, key)
            @warn "CharmmFF: ignoring unknown force field option $key"
            continue
        end
        options[key] = options[key] isa T ? T(value) : value
    end

    charmm_ff = ForceField{T}(
        "CharmmFF",
        ac,
        charmm_params,
        options,
        init_atom_types(charmm_params, T),
        Vector{AbstractForceFieldComponent{T}}(),
        Dict{String, T}(),
        Vector{Atom{T}}(),
        constrained_atoms
    )

    if options[:remove_extended_hydrogens]
        _remove_extended_hydrogens!(charmm_ff)
    end

    assign_typenames_and_charges!(charmm_ff)

    append!(
        charmm_ff.components,
        [
            QuadraticStretchComponent{T}(charmm_ff),
            QuadraticBendComponent{T}(charmm_ff),
            CharmmTorsionComponent{T}(charmm_ff),
            CharmmImproperTorsionComponent{T}(charmm_ff),
            CharmmNonBondedComponent{T}(charmm_ff)
        ]
    )

    setup!(charmm_ff)
    update!(charmm_ff)
    _check_warnings(charmm_ff)

    charmm_ff
end
