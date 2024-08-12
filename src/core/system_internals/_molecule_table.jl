const _molecule_table_schema = Tables.Schema(
    (:idx, :name),
    (Int, String)
)
const _molecule_table_cols = _molecule_table_schema.names
const _molecule_table_cols_set = Set(_molecule_table_cols)
const _molecule_table_cols_priv = Set([:variant, :properties, :flags])

@auto_hash_equals struct _MoleculeTable <: AbstractColumnTable
    # public columns
    idx::Vector{Int}
    name::Vector{String}

    # private columns
    variant::Vector{MoleculeVariantType}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # internals
    _idx_map::Dict{Int,Int}

    function _MoleculeTable()
        new(
            Int[],
            String[],
            MoleculeVariantType[],
            Properties[],
            Flags[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.columnnames(::_MoleculeTable) = _molecule_table_cols
@inline Tables.schema(::_MoleculeTable) = _molecule_table_schema

@inline function Tables.getcolumn(mt::_MoleculeTable, nm::Symbol)
    @assert nm in _molecule_table_cols_priv || nm in _molecule_table_cols "type _MoleculeTable has no column $nm"
    getfield(mt, nm)
end

@inline Base.size(mt::_MoleculeTable) = (length(mt.idx), length(_molecule_table_cols))

function Base.push!(
    mt::_MoleculeTable,
    idx::Int = 0;
    name::String = "",
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

@inline _row_by_idx(mt::_MoleculeTable, idx::Int) = ColumnTableRow(getfield(mt, :_idx_map)[idx], mt)
