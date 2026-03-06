export
    FragmentDB,
    default_fragmentdb,
    infer_topology!

@auto_hash_equals struct DBAtom{T <: Real}
    name::String
    element::ElementType
    r::Vector{T}
end

@auto_hash_equals struct DBBond
    number::Int
    a1::String
    a2::String
    order::BondOrderType
end

@auto_hash_equals struct DBConnection{T <: Real}
    name::String
    atom_name::String
    match_name::String
    order::BondOrderType
    distance::T
    tolerance::T
end

@auto_hash_equals struct DBVariant
    name::String
    delete::Union{Nothing, Set{String}}
    rename::Union{Nothing, Dict{String, String}}
    flags::Union{Nothing, Flags}
end

@auto_hash_equals struct DBFragment{T <: Real}
    name::String
    names::Vector{String}
    type::FragmentVariantType
    atoms::Vector{DBAtom{T}}
    bonds::Vector{DBBond}
    connections::Union{Nothing, Vector{DBConnection{T}}}
    variants::Union{Nothing, Vector{DBVariant}}
end

@auto_hash_equals struct DBFragmentVariant{T <: Real}
    name::String
    type::FragmentVariantType
    atoms::Vector{DBAtom{T}}
    bonds::Vector{DBBond}
    flags::Union{Nothing, Flags}

    function DBFragmentVariant{T}(frag::DBFragment{T}, var::DBVariant) where T
        atoms = frag.atoms
        bonds = frag.bonds

        if !isnothing(var.delete)
            atoms = filter(a -> a.name ∉ var.delete, atoms)
            bonds = filter(b -> b.a1 ∉ var.delete && b.a2 ∉ var.delete, bonds)
        end

        if !isnothing(var.rename)
            atoms = map(a -> DBAtom{T}(get(var.rename, a.name, a.name), a.element, a.r), atoms)
            bonds = map(b -> DBBond(b.number, get(var.rename, b.a1, b.a1), get(var.rename, b.a2, b.a2), b.order), bonds)
        end

        new{T}(var.name, frag.type, atoms, bonds, var.flags)
    end

    @inline function DBFragmentVariant{T}(frag::DBFragment{T}) where T
        new{T}(frag.name, frag.type, frag.atoms, frag.bonds, nothing)
    end
end

@auto_hash_equals struct DBNameMapping
    name::String
    maps_to::String
    mappings::Dict{String, String}
end

"""
    $(TYPEDEF)

Fragment database. Contains information about common biomolecule fragments such as amino acids
and nucleotides that can be used to reconstruct incomplete systems.

# Constructors
    FragmentDB(path::AbstractString = ball_data_path("fragmentdb/FragmentDB.json"))

Creates a new `FragmentDB{Float32}` from the given configuration file.

    FragmentDB{T}(path::AbstractString = ball_data_path("fragmentdb/FragmentDB.json"))

Creates a new `FragmentDB{T}` from the given configuration file.
"""
@auto_hash_equals struct FragmentDB{T <: Real}
    fragments::OrderedDict{String, DBFragment{T}}
    name_mappings::OrderedDict{String, DBNameMapping}
    defaults::OrderedDict{String, String}

    function FragmentDB{T}(path::AbstractString = ball_data_path("fragmentdb/FragmentDB.json")) where T
        dir = dirname(path)
        fdb = JSON.parse(read(path, String))
        fragments = OrderedDict(
            e.first => JSON.parse(read(joinpath(dir, e.second), String), DBFragment{T})
            for e in fdb.fragments
        )
        name_mappings = OrderedDict(
            e.first => JSON.parse(read(joinpath(dir, e.second), String), DBNameMapping)
            for e in fdb.names
        )
        defaults = OrderedDict(e.first => e.second for e in fdb.defaults)
        new{T}(fragments, name_mappings, defaults)
    end
end

@inline function FragmentDB(path::AbstractString = ball_data_path("fragmentdb/FragmentDB.json"))
    FragmentDB{Float32}(path)
end

"""
    const _default_fragmentdb

Global default fragment database.
"""
const _default_fragmentdb = FragmentDB()

"""
    $(TYPEDSIGNATURES)

Returns the global default fragment database.
"""
@inline function default_fragmentdb()
    _default_fragmentdb
end

@inline function Base.show(io::IO, fdb::FragmentDB)
    print(io,
        "FragmentDB with ",
        length(fdb.fragments), " fragments, ",
        length(fdb.name_mappings), " mappings, ",
        length(fdb.defaults), " defaults."
    )
end

@inline function bonds(a::DBAtom{T}, var::DBFragmentVariant{T}) where T
    filter(b -> a.name ∈ [b.a1, b.a2], var.bonds)
end

@inline function get_partner(bond::DBBond, atom::DBAtom{T}, var::DBFragmentVariant{T}) where T
    if bond.a1 == atom.name
        return var.atoms[findfirst(a -> a.name == bond.a2, var.atoms)]
    elseif bond.a2 == atom.name
        return var.atoms[findfirst(a -> a.name == bond.a1, var.atoms)]
    else
        return nothing
    end
end

function label_terminal_fragments!(ac::AbstractAtomContainer{T}, fdb::FragmentDB{T}) where T
    # iterate over all chains and label their terminals
    for chain in chains(ac)
        ft = fragments(chain)
        length(ft) == 0 && continue

        # fragment variants are possibly not assigned yet, so we rely
        # on the FragmentDB
        is_amino_acid = f -> haskey(fdb.fragments, f.name) &&
            fdb.fragments[f.name].type == FragmentVariant.Residue
        is_nucleotide = f -> haskey(fdb.fragments, f.name) &&
            fdb.fragments[f.name].type == FragmentVariant.Nucleotide

        # first, the n- and c-terminals
        amino_acids = filter(is_amino_acid, ft)
        if length(amino_acids) > 0
            set_flag!(first(amino_acids), :N_TERMINAL)
            set_flag!(last(amino_acids),  :C_TERMINAL)
        end

        # then, the 5'- and 3'-primes (nucleic acid chains run 5' → 3')
        nucleotides = filter(is_nucleotide, ft)
        if length(nucleotides) > 0
            set_flag!(first(nucleotides), Symbol("5_PRIME"))
            set_flag!(last(nucleotides),  Symbol("3_PRIME"))
        end
    end
    nothing
end

@inline function label_terminal_fragments!(ac::AbstractAtomContainer{Float32})
    label_terminal_fragments!(ac, default_fragmentdb())
end

function get_reference_fragment(f::Fragment{T}, fdb::FragmentDB{T}) where T
    # first, try to find the fragment in the database
    if f.name ∉ keys(fdb.fragments)
        return nothing
    end

    db_fragment = fdb.fragments[f.name]

    # fragments w/o declared variants, e.g., skeletons
    if isnothing(db_fragment.variants)
        return DBFragmentVariant{T}(db_fragment)
    end

    # does the fragment have variants?
    if length(db_fragment.variants) == 1
        return DBFragmentVariant{T}(db_fragment, only(db_fragment.variants))
    end

    # now, find the variant that best matches the fragment
    # This returns N/C terminal variants for fragments that
    # have the corresponding flags set or cysteine variants
    # without thiol hydrogen if the `HAS_SSBOND` flag is set
    best_score = -1
    best_variant = nothing

    # iterate over all variants of the fragment and compare the flags
    for var in db_fragment.variants
        vflags = isnothing(var.flags) ? Flags() : var.flags

        # count the flags fragment and variant have in common
        score = length(f.flags ∩ vflags)

        # subtract the number of variant props that the fragment doesn't have
        score -= length(setdiff(vflags, f.flags))

        @debug "Considering variant $(var.name). # score: $score"

        if (score > best_score)
            best_variant = var
            best_score = score
        end
    end

    if isnothing(best_variant)
        return nothing
    end

    DBFragmentVariant{T}(db_fragment, best_variant)
end

@inline function get_reference_fragment(f::Fragment{Float32})
    get_reference_fragment(f, default_fragmentdb())
end

"""
    infer_topology!(::AbstractAtomContainer{Float32})
    infer_topology!(::AbstractAtomContainer{T}, ::FragmentDB{T})

Apply standard preprocessing functions to the given atom container, using the default/given
fragment database. By default, this calls (in order): [`normalize_names!`](@ref),
[`reconstruct_fragments!`](@ref), and [`build_bonds!`](@ref).

# Supported keyword arguments
 - `normalize_names::Bool = true`
 - `reconstruct_fragments::Bool = true`
 - `build_bonds::Bool = true`
All keyword arguments enable or disable the corresponding preprocessors.

!!! warning
    Caution is advised when partly enabling/disabling preprocessors, as they might depend on
    each other. Please refer to the corresponding documentation.
"""
@inline function infer_topology!(
    ac::AbstractAtomContainer{T},
    fdb::FragmentDB{T};
    normalize_names::Bool = true,
    reconstruct_fragments::Bool = true,
    build_bonds::Bool = true
) where T
    normalize_names       && normalize_names!(ac, fdb)
    reconstruct_fragments && reconstruct_fragments!(ac, fdb)
    build_bonds           && build_bonds!(ac, fdb)
    nothing
end

@inline function infer_topology!(ac::AbstractAtomContainer{Float32})
    infer_topology!(ac, default_fragmentdb())
end
