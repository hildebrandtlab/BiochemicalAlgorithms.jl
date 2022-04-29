using NamedTupleTools

const ProteinAtomExtras = @NamedTuple begin
    residue_id::Int
    residue_name::String
    chain::String
end

const ProteinAtom{T} = merge(Atom{T}, ProteinAtomExtras)