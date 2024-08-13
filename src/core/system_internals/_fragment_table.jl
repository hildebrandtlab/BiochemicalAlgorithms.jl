const _fragment_table_schema = Tables.Schema(
    (:idx, :number, :name),
    (Int, Int, String)
)
const _fragment_table_cols = _fragment_table_schema.names
const _fragment_table_cols_set = Set(_fragment_table_cols)
const _fragment_table_cols_priv = Set([:variant, :properties, :flags, :molecule_idx, :chain_idx])

@auto_hash_equals struct _FragmentTable <: AbstractColumnTable
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    name::Vector{String}

    # private columns
    variant::Vector{FragmentVariantType}
    properties::Vector{Properties}
    flags::Vector{Flags}
    molecule_idx::Vector{Int}
    chain_idx::Vector{Int}

    # internals
    _idx_map::Dict{Int,Int}

    function _FragmentTable()
        new(
            Int[],
            Int[],
            String[],
            FragmentVariantType[],
            Properties[],
            Flags[],
            Int[],
            Int[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.columnnames(::_FragmentTable) = _fragment_table_cols
@inline Tables.schema(::_FragmentTable) = _fragment_table_schema

@inline function Tables.getcolumn(ft::_FragmentTable, nm::Symbol)
    @assert nm in _fragment_table_cols_priv || nm in _fragment_table_cols_set "type _FragmentTable has no column $nm"
    getfield(ft, nm)
end

@inline Base.size(ft::_FragmentTable) = (length(ft.idx), length(_fragment_table_cols))

function Base.push!(
    ft::_FragmentTable,
    idx::Int,
    number::Int,
    molecule_idx::Int,
    chain_idx::Int;
    name::String = "",
    variant::FragmentVariantType = FragmentVariant.None,
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    ft._idx_map[idx] = length(ft.idx) + 1
    push!(ft.idx, idx)
    push!(ft.number, number)
    push!(ft.name, name)
    push!(ft.variant, variant)
    push!(ft.properties, properties)
    push!(ft.flags, flags)
    push!(ft.molecule_idx, molecule_idx)
    push!(ft.chain_idx, chain_idx)
    ft
end

function _fragment_table(itr)
    ft = _FragmentTable()
    for f in itr
        push!(ft, f.idx, f.number, f.molecule_idx, f.chain_idx;
            name = f.name,
            variant = f.variant,
            properties = f.properties,
            flags = f.flags
       )
    end
    ft
end
@inline Tables.materializer(::Type{_FragmentTable}) = itr -> _fragment_table(itr)

@inline _rowno_by_idx(ft::_FragmentTable, idx::Int) = getindex(getfield(ft, :_idx_map), idx)
@inline _row_by_idx(ft::_FragmentTable, idx::Int) = ColumnTableRow(_rowno_by_idx(ft, idx), ft)
