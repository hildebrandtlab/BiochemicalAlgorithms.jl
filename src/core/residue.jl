export Residue

const Residue = @NamedTuple begin
    number::Int
    type::AminoAcid
    chain_id::Int
end