export
    AmberFF

get_amber_default_options(T=Float32) = Dict{Symbol, Any}(
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

function AmberFF(
        ac::AbstractAtomContainer{T},
        filename=ball_data_path("forcefields/AMBER/amber96.ini");
        constrained_atoms=Vector{Int}()) where {T<:Real}

    amber_params = AmberFFParameters(filename)
    amber_ff = ForceField{T}(
        "AmberFF",
        ac,
        amber_params,
        get_amber_default_options(T),
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

    amber_ff
end
