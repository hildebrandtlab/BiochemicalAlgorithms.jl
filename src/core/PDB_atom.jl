using NamedTupleTools

const PDBAtomExtras = @NamedTuple begin
    residue_id::Int
    residue_name::String
    chain::String
end

const PDBAtom{T} = merge(Atom{T}, PDBAtomExtras)