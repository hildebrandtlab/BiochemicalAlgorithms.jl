export
    AmberFF

get_amber_default_options(::Type{T} = Float32) where {T <: Real} = Dict{Symbol, Any}(
    :nonbonded_cutoff               => T(20.0),
    :vdw_cutoff                     => T(15.0),
    :vdw_cuton                      => T(13.0),
    :electrostatic_cutoff           => T(15.0),
    :electrostatic_cuton            => T(13.0),
    :scaling_vdw_1_4                => T(2.0),
    :scaling_electrostatic_1_4      => T(1.2),
    :distance_dependent_dielectric  => false,
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

"""
    AmberFF(
        ::AbstractAtomContainer{T},
        param_file::AbstractString = ball_data_path("forcefields/AMBER/amber96.ini")
    )

Initializes an AMBER force field for the given atom container and the given parameter
file (default: AMBER96).

# Supported keyword arguments
 - `nonbonded_cutoff::T = 20`
 - `vdw_cutoff::T = 15`
 - `vdw_cuton::T = 13`
 - `electrostatic_cutoff::T = 15`
 - `electrostatic_cuton::T = 13`
 - `scaling_vdw_1_4 = 2`
 - `scaling_electrostatic_1_4::T = 1.2`
 - `distance_dependent_dielectric::Bool = false`
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
function AmberFF(
    ac::AbstractAtomContainer{T},
    param_file::AbstractString=ball_data_path("forcefields/AMBER/amber96.ini");
    constrained_atoms=Vector{Int}(),
    kwargs...
) where T

    amber_params = AmberFFParameters(param_file, T)
    options = get_amber_default_options(T)

    # read :scaling_electrostatic_1_4 from parameter file (if present)
    if haskey(amber_params.sections, "Options")
        props = extract_section(amber_params, "Options").properties
        if haskey(props, "SCEE")
            options[:scaling_electrostatic_1_4] = parse(T, props["SCEE"])
        end
    end

    for (key, value) in kwargs
        if !haskey(options, key)
            @warn "AmberFF: ignoring unknown force field option $key"
            continue
        end
        options[key] = value
    end

    amber_ff = ForceField{T}(
        "AmberFF",
        ac,
        amber_params,
        options,
        init_atom_types(amber_params, T),
        Vector{AbstractForceFieldComponent{T}}(),
        Dict{String, T}(),
        Vector{Atom{T}}(),
        constrained_atoms
    )

    assign_typenames_and_charges!(amber_ff)

    append!(
        amber_ff.components,
        [
            QuadraticStretchComponent{T}(amber_ff),
            QuadraticBendComponent{T}(amber_ff),
            TorsionComponent{T}(amber_ff),
            NonBondedComponent{T}(amber_ff)
        ]
    )

    setup!(amber_ff)
    update!(amber_ff)
    _check_warnings(amber_ff)

    amber_ff
end
