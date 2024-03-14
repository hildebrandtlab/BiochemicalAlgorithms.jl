export
    AbstractColumnTable,
    AbstractSystemComponent,
    AbstractSystemComponentTable,
    get_property,
    has_flag,
    has_property,
    set_flag!,
    set_property!,
    unset_flag!

"""
    $(TYPEDEF)

Abstract base type for all Tables.jl-compatible column tables.
"""
abstract type AbstractColumnTable <: Tables.AbstractColumns end

@inline Tables.istable(::Type{<: AbstractColumnTable}) = true
@inline Tables.columnaccess(::Type{<: AbstractColumnTable}) = true
@inline Tables.columns(at::AbstractColumnTable) = at
@inline Tables.rows(at::AbstractColumnTable) = at
@inline Tables.getcolumn(at::AbstractColumnTable, i::Int) = Tables.getcolumn(at, Tables.columnnames(at)[i])

@inline Base.getproperty(at::AbstractColumnTable, nm::Symbol) = getfield(at, nm)
@inline Base.size(at::AbstractColumnTable, dim) = size(at)[dim]
@inline Base.length(at::AbstractColumnTable) = size(at, 1)
@inline Base.keys(at::AbstractColumnTable) = LinearIndices((length(at),))

@inline function Base.show(io::IO, at::AbstractColumnTable) 
    print(io, "$(typeof(at)) with $(length(at)) rows")
end

@inline function Base.show(io::IO, ::MIME"text/html", at::AbstractColumnTable)
    show(io, at; backend = Val(:html))
end

@inline function Base.show(io::IO, ::MIME"text/plain", at::AbstractColumnTable)
    _show(io, at; reserved_display_lines = 2)
end

@inline function _show(io::IO, at::AbstractColumnTable; kwargs...)
    PrettyTables.pretty_table(io, at;
        alignment = :l,
        header = collect(Tables.columnnames(at)),
        max_num_of_rows = 50,
        row_number_column_title = "#",
        show_row_number = true,
        title  = "$(typeof(at)) with $(length(at)) rows:",
        vcrop_mode = :middle,
        kwargs...
    )
end

@inline function _filter_select(itr, col::Symbol)
    (getproperty(a, col) for a in itr)
end

"""
    $(TYPEDEF)

Abstract base type for all Tables.jl-compatible system component tables.
"""
abstract type AbstractSystemComponentTable{T <: Real} <: AbstractColumnTable end

"""
    $(TYPEDEF)

Abstract base type for all components of a system, including the system itself.
"""
abstract type AbstractSystemComponent{T <: Real} end

#=
    Properties
=#

"""
    $(TYPEDSIGNATURES)

Returns a `Bool` indicating whether the given system component has the given property.
"""
@inline has_property(ac::AbstractSystemComponent, key::Symbol) = haskey(ac.properties, key)

"""
    $(TYPEDSIGNATURES)

Returns the property associated with the given key in `ac`.
"""
@inline get_property(ac::AbstractSystemComponent, key::Symbol) = ac.properties[key]

"""
    $(TYPEDSIGNATURES)

Returns the property associated with the given key in `ac`. If no such property exists, returns `default`.
"""
@inline get_property(ac::AbstractSystemComponent, key::Symbol, default) = get(ac.properties, key, default)

"""
    $(TYPEDSIGNATURES)

Sets the property associated with the given key in `ac` to the given `value`.
"""
@inline set_property!(ac::AbstractSystemComponent, key::Symbol, value) = ac.properties[key] = value

#=
    Flags
=#

"""
    $(TYPEDSIGNATURES)

Returns a `Bool` indicating whether the given system component has the given flag.
"""
@inline has_flag(ac::AbstractSystemComponent, flag::Symbol) = flag in ac.flags

"""
    $(TYPEDSIGNATURES)

Adds the given flag to `ac`.
"""
@inline set_flag!(ac::AbstractSystemComponent, flag::Symbol) = push!(ac.flags, flag)

"""
    $(TYPEDSIGNATURES)

Removes the given flag from `ac`.
"""
@inline unset_flag!(ac::AbstractSystemComponent, flag::Symbol) = delete!(ac.flags, flag)
