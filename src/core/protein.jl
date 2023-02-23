using DataFrames

export Protein

const _ProteinAtom{T} = NamedTuple{
    (fieldnames(Atom{T})..., :frame_id, :residue_id),
    Tuple{fieldtypes(Atom{T})..., Int, Int}
}

mutable struct Protein{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame
    residues::DataFrame
    properties::Properties

    function Protein{T}(name = "") where T 
        new(name, DataFrame(_ProteinAtom{T}[]), DataFrame(Bond[]), DataFrame(Residue[]), Properties())
    end
end

Protein(name="") = Protein{Float32}(name)

@inline function residues(ac::Protein{T}) where T
    ac.residues
end

@inline function atoms(ac::Protein{T}; frame_id::Int = 1) where T
    df = get(groupby(ac.atoms, :frame_id), (frame_id = frame_id,), DataFrame(_ProteinAtom{T}[]))
    view(df, :, 1:length(fieldnames(Atom{T})))
end

@inline function atoms(ac::Protein{T}, r::Residue; frame_id::Int = 1) where T
    df = get(
        groupby(ac.atoms, [:frame_id, :residue_id]),
        (frame_id = frame_id, residue_id = r.number),
        DataFrame(_ProteinAtom{T}[])
    )
    view(df, :, 1:length(fieldnames(Atom{T})))
end

function Base.push!(ac::Protein{T}, r::Residue, atoms::Atom{T}...; frame_id::Int = 1) where T
    # check whether the given residue is known, add it otherwise
    res_df = ac.residues
    if isempty(@view res_df[res_df.number .== r.number, :])
        push!(res_df, r)
    end

    # add all atoms
    for atom in atoms
        push!(ac.atoms, (atom..., frame_id = frame_id, residue_id = r.number))
    end
end
