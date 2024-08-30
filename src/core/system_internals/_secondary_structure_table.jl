const _secondary_structure_table_cols_main = (:idx, :number, :type, :name)
const _secondary_structure_table_cols_extra = (:properties, :flags, :molecule_idx, :chain_idx)
const _secondary_structure_table_cols = (_secondary_structure_table_cols_main..., _secondary_structure_table_cols_extra...)
const _secondary_structure_table_cols_set = Set(_secondary_structure_table_cols)
const _secondary_structure_table_schema = Tables.Schema(
    _secondary_structure_table_cols_main,
    (Int, Int, SecondaryStructureType, String)
)

@auto_hash_equals struct _SecondaryStructureTable <: _AbstractSystemComponentTable
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    type::Vector{SecondaryStructureType}
    name::Vector{String}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}
    molecule_idx::Vector{Int}
    chain_idx::Vector{Int}

    # internals
    _idx_map::_IdxMap

    function _SecondaryStructureTable()
        new(
            Int[],
            Int[],
            SecondaryStructureType[],
            String[],
            Properties[],
            Flags[],
            Int[],
            Int[],
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_SecondaryStructureTable, Type{_SecondaryStructureTable}}, nm::Symbol)
    nm in _secondary_structure_table_cols_set
end

@inline function Tables.columnnames(::_SecondaryStructureTable)
    _secondary_structure_table_cols_main
end

@inline function Tables.schema(::_SecondaryStructureTable)
    _secondary_structure_table_schema
end

@inline function Base.propertynames(::_SecondaryStructureTable)
    _secondary_structure_table_cols
end

function Base.push!(
    st::_SecondaryStructureTable,
    idx::Int,
    number::Int,
    type::SecondaryStructureType,
    molecule_idx::Int,
    chain_idx::Int;
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    st._idx_map[idx] = length(st.idx) + 1
    push!(st.idx, idx)
    push!(st.number, number)
    push!(st.type, type)
    push!(st.name, name)
    push!(st.properties, properties)
    push!(st.flags, flags)
    push!(st.molecule_idx, molecule_idx)
    push!(st.chain_idx, chain_idx)
    st
end

function _secondary_structure_table(itr)
    st = _SecondaryStructureTable()
    for c in itr
        push!(st, c.idx, c.number, c.type, c.molecule_idx, c.chain_idx;
            name = c.name,
            properties = c.properties,
            flags = c.flags
       )
    end
    st
end
@inline Tables.materializer(::Type{_SecondaryStructureTable}) = itr -> _secondary_structure_table(itr)
