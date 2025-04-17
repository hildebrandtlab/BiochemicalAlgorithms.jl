@auto_hash_equals struct AminoAcidDescriptor{T<:Real}
    type::String
    angle_phi::T
    angle_psi::T
    angle_omega::T

    function AminoAcidDescriptor(type::String, angle_phi::T, angle_psi::T, angle_omega::T) where {T<:Real}
        new{T}(type, angle_phi, angle_psi, angle_omega)
    end
end

function PeptideBuilder(amino_acid_descriptors::Vector{AminoAcidDescriptor{T}}) where {T<:Real}
    new{T}(amino_acid_descriptors)
    amino_acid_descriptors::Vector{AminoAcidDescriptor{T}}
end

#function addAminoAcid!(builder::PeptideBuilder{T}, type::String, angle_phi::T, angle_psi::T, angle_omega::T) where {T<:Real}
#    push!(builder.amino_acid_descriptors, AminoAcidDescriptor(type, angle_phi, angle_psi, angle_omega))
#end