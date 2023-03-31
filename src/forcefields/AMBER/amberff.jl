export AmberFF

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
    :max_number_of_unassigned_atoms => typemax(Int32)
)

function AmberFF(
        ac::AbstractAtomContainer{T}, 
        filename=ball_data_path("forcefields/AMBER/amber96.ini")) where {T<:Real}

    amber_params = AmberFFParameters(filename)
    amber_ff = ForceField{T}(
        "AmberFF",
        ac, 
        amber_params, 
        get_amber_default_options(T),
        init_atom_types(amber_params),
        Vector{AbstractForceFieldComponent{T}}(),
        Vector{T}(),
        Vector{Atom{T}}()
    )

    assign_typenames_and_charges!(amber_ff)

    push!(
        amber_ff.components,
        QuadraticStretchComponent{T}(amber_ff)
    )

    push!(
        amber_ff.energies,
        zero(T)
    )

    amber_ff
end