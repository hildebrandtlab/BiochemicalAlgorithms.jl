const _molecule_table_cols_main = (:idx, :name)
const _molecule_table_cols_extra = (:variant, :properties, :flags)
const _molecule_table_cols = (_molecule_table_cols_main..., _molecule_table_cols_extra...)
const _molecule_table_cols_set = Set(_molecule_table_cols)
const _molecule_table_schema = Tables.Schema(
    _molecule_table_cols_main,
    (Int, String)
)

@auto_hash_equals struct _MoleculeTable <: _AbstractSystemComponentTable
    # public columns
    idx::Vector{Int}
    name::Vector{String}

    # private columns
    variant::Vector{MoleculeVariantType}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # internals
    _idx_map::_IdxMap

    function _MoleculeTable()
        new(
            Int[],
            String[],
            MoleculeVariantType[],
            Properties[],
            Flags[],
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_MoleculeTable, Type{_MoleculeTable}}, nm::Symbol)
    nm in _molecule_table_cols_set
end

@inline function Tables.columnnames(::_MoleculeTable)
    _molecule_table_cols_main
end

@inline function Tables.schema(::_MoleculeTable)
    _molecule_table_schema
end

@inline function Base.propertynames(::_MoleculeTable)
    _molecule_table_cols
end

function Base.push!(
    mt::_MoleculeTable,
    idx::Int = 0;
    name::AbstractString = "",
    variant::MoleculeVariantType = MoleculeVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    mt._idx_map[idx] = length(mt.idx) + 1
    push!(mt.idx, idx)
    push!(mt.name, name)
    push!(mt.variant, variant)
    push!(mt.properties, properties)
    push!(mt.flags, flags)
    mt
end

function _molecule_table(itr)
    mt = _MoleculeTable()
    for m in itr
        push!(mt, m.idx;
            name = m.name,
            variant = m.variant,
            properties = m.properties,
            flags = m.flags
        )
    end
    mt
end
Tables.materializer(::Type{_MoleculeTable}) = _molecule_table
