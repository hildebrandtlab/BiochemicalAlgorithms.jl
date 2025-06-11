export
    AminoAcidDescriptor,
    PeptideBuilder

"""
    Immutable representation of an individual amino acid in the sequence.

# Public fields
- `type::String`
- `angle_phi::T`
- `angle_psi::T`
- `angle_omega::T`

# Private fields


# Constructors


"""
@auto_hash_equals struct AminoAcidDescriptor{T}
    type::String
    angle_phi::T
    angle_psi::T
    angle_omega::T

    function AminoAcidDescriptor(
        type::String,
        angle_phi::T=deg2rad(T(-47)),
        angle_psi::T=deg2rad(T(-58)),
        angle_omega::T=deg2rad(T(180))
    ) where {T}
        new{T}(type, angle_phi, angle_psi, angle_omega) #TODO: check if type is a valid amino acid, switch to 3letter code
    end
end

@auto_hash_equals struct PeptideBuilder{T}
    sequence::Vector{AminoAcidDescriptor{T}}
    chain_name::String
    protein_name::String
    is_proline::Bool
    fragment_db::FragmentDB

    function PeptideBuilder(
        sequence::Vector{AminoAcidDescriptor{T}}=AminoAcidDescriptor(),
        chain_name::String="A",
        protein_name::String="",
        is_proline::Bool=false,
        fragment_db::FragmentDB=FragmentDB()
    ) where {T}
        new{T}(sequence, chain_name, protein_name, is_proline, fragment_db)
    end

end

function addAminoAcid!(builder::PeptideBuilder{T}, type::String, angle_phi::T, angle_psi::T, angle_omega::T) where {T<:Real}
    push!(builder.amino_acid_descriptors, AminoAcidDescriptor(type, angle_phi, angle_psi, angle_omega))
end

