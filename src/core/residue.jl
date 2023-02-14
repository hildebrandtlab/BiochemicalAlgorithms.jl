export Residue

const Residue = @NamedTuple begin
    number::Int
    type::AminoAcid
    chain::String
end
