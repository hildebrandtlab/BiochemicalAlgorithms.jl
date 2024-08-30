const _fragment_table_cols_main = (:idx, :number, :name)
const _fragment_table_cols_extra = (:variant, :properties, :flags, :molecule_idx, :chain_idx, :secondary_structure_idx)
const _fragment_table_cols = (_fragment_table_cols_main..., _fragment_table_cols_extra...)
const _fragment_table_cols_set = Set(_fragment_table_cols)
const _fragment_table_schema = Tables.Schema(
    _fragment_table_cols_main,
    (Int, Int, String)
)

@auto_hash_equals struct _FragmentTable <: _AbstractSystemComponentTable
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
    secondary_structure_idx::Vector{MaybeInt}

    # internals
    _idx_map::_IdxMap

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
            MaybeInt[],
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_FragmentTable, Type{_FragmentTable}}, nm::Symbol)
    nm in _fragment_table_cols_set
end

@inline function Tables.columnnames(::_FragmentTable)
    _fragment_table_cols_main
end

@inline function Tables.schema(::_FragmentTable)
    _fragment_table_schema
end

@inline function Base.propertynames(::_FragmentTable)
    _fragment_table_cols
end

function Base.push!(
    ft::_FragmentTable,
    idx::Int,
    number::Int,
    molecule_idx::Int,
    chain_idx::Int;
    secondary_structure_idx::MaybeInt = nothing,
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
    push!(ft.secondary_structure_idx, secondary_structure_idx)
    ft
end

function _fragment_table(itr)
    ft = _FragmentTable()
    for f in itr
        push!(ft, f.idx, f.number, f.molecule_idx, f.chain_idx;
            secondary_structure_idx = f.secondary_structure_idx,
            name = f.name,
            variant = f.variant,
            properties = f.properties,
            flags = f.flags
       )
    end
    ft
end
@inline Tables.materializer(::Type{_FragmentTable}) = itr -> _fragment_table(itr)
