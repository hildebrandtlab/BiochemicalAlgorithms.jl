export PDBAtom

const PDBAtomExtras = @NamedTuple begin
    residue_id::Int
    residue_name::String
    chain::String
end

const PDBAtom{T} = NamedTuple{
    (fieldnames(Atom{T})..., fieldnames(PDBAtomExtras)...), 
    Tuple{fieldtypes(Atom{T})..., fieldtypes(PDBAtomExtras)...}
}
