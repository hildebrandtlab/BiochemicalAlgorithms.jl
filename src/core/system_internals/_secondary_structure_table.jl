const _secondary_structure_table_schema = Tables.Schema(
    (:idx, :number, :type, :name),
    (Int, Int, SecondaryStructureType, String)
)
const _secondary_structure_table_cols = _secondary_structure_table_schema.names
const _secondary_structure_table_cols_set = Set(_secondary_structure_table_cols)
const _secondary_structure_table_cols_priv = Set([:properties, :flags, :molecule_idx, :chain_idx])

@auto_hash_equals struct _SecondaryStructureTable <: AbstractColumnTable
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
    _idx_map::Dict{Int,Int}

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
            Dict{Int,Int}()
        )
    end
end

@inline Tables.columnnames(::_SecondaryStructureTable) = _secondary_structure_table_cols
@inline Tables.schema(::_SecondaryStructureTable) = _secondary_structure_table_schema

@inline function Tables.getcolumn(st::_SecondaryStructureTable, nm::Symbol)
    @assert nm in _secondary_structure_table_cols_priv || nm in _secondary_structure_table_cols_set "type _SecondaryStructureTable has no column $nm"
    getfield(st, nm)
end

@inline Base.size(st::_SecondaryStructureTable) = (length(st.idx), length(_secondary_structure_table_cols))

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

@inline function _rebuild_idx_map!(st::_SecondaryStructureTable)
    empty!(st._idx_map)
    merge!(st._idx_map, Dict(v => k for (k, v) in enumerate(st.idx)))
    st
end

function _delete!(st::_SecondaryStructureTable, rowno::Int)
    deleteat!(st.idx, rowno)
    deleteat!(st.number, rowno)
    deleteat!(st.type, rowno)
    deleteat!(st.name, rowno)
    deleteat!(st.properties, rowno)
    deleteat!(st.flags, rowno)
    deleteat!(st.molecule_idx, rowno)
    deleteat!(st.chain_idx, rowno)
    nothing
end

function Base.delete!(st::_SecondaryStructureTable, idx::Int)
    _delete!(st, st._idx_map[idx])
    _rebuild_idx_map!(st)
end

function Base.delete!(st::_SecondaryStructureTable, idx::Vector{Int})
    rownos = getindex.(Ref(st._idx_map), idx)
    unique!(rownos)
    sort!(rownos; rev = true)
    for rowno in rownos
        _delete!(st, rowno)
    end
    _rebuild_idx_map!(st)
end

function Base.empty!(st::_SecondaryStructureTable)
    empty!(st.idx)
    empty!(st.number)
    empty!(st.type)
    empty!(st.name)
    empty!(st.properties)
    empty!(st.flags)
    empty!(st.molecule_idx)
    empty!(st.chain_idx)
    empty!(st._idx_map)
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

@inline _rowno_by_idx(st::_SecondaryStructureTable, idx::Int) = getindex(getfield(st, :_idx_map), idx)
@inline _row_by_idx(st::_SecondaryStructureTable, idx::Int) = ColumnTableRow(_rowno_by_idx(st, idx), st)
