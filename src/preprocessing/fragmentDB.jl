using StructTypes
using JSON3
using AutoHashEquals

import DataStructures: OrderedCollections.OrderedDict

export FragmentDB

@auto_hash_equals struct DBNode
    key::String
    value::Union{AbstractArray{DBNode}, String}
end
StructTypes.StructType(::Type{DBNode}) = StructTypes.Struct()

@auto_hash_equals struct DBAtom{T<:Real}
    name::String
    element::ElementType
    r::Vector3

    function DBAtom{T}(n::DBNode) where {T<:Real}
        name = n.key

        raw_data = split(n.value)

        if length(raw_data) == 4
            element = parse(Elements, raw_data[1])

            r = Vector3(map(d -> parse(T, d), raw_data[2:4]))

            return new(name, element, r)
        end

        throw(ArgumentError("DBAtom: invalid format!"))
    end
end
StructTypes.StructType(::Type{DBAtom{T}}) where {T<:Real} = StructTypes.Struct()

DBAtom(n::DBNode) = DBAtom{Float32}(n)

function find_atom(name, atoms::Vector{DBAtom{T}}) where {T<:Real}
    candidates = filter(a -> a.name == name, atoms)

    if length(candidates) != 1
        throw(ArgumentError("FragmentDB::find_atom: atom not found!"))
    end
    
    candidates[1]
end

function to_bond_order(s)
    order = BondOrder.Unknown

    if s == "s"
        order = BondOrder.Single
    elseif s == "d"
        order = BondOrder.Double
    elseif s == "t"
        order = BondOrder.Triple
    elseif s == "a"
        order = BondOrder.Aromatic
    else
        throw(ArgumentError("DBBond: invalid format!"))
    end

    order
end

@auto_hash_equals struct DBBond{T<:Real}
    number::Int
    a1::String
    a2::String
    order::BondOrderType

    function DBBond{T}(n::DBNode) where {T<:Real}
        number = parse(Int, n.key)

        raw_data = split(n.value)

        a1 = raw_data[1]
        a2 = raw_data[2]
        
        order = to_bond_order(raw_data[3])

        new(number, a1, a2, order)
    end
end

DBBond(n::DBNode) = DBond{Float32}(n)

@auto_hash_equals struct DBConnection{T<:Real}
    name::String
    atom_name::String
    match_name::String
    order::BondOrderType
    distance::T
    tolerance::T

    function DBConnection{T}(n::DBNode) where {T<:Real}
        name = n.key

        raw_data = split(n.value)

        if length(raw_data) != 5
            throw(ArgumentError("DBConnection: invalid format!"))
        end

        atom_name  = raw_data[1]
        match_name = raw_data[2]
        order      = to_bond_order(raw_data[3])

        distance   = parse(T, raw_data[4])
        tolerance  = parse(T, raw_data[5])

        new(name, atom_name, match_name, order, distance, tolerance)
    end
end

DBConnection(n::DBNode) = DBConnection{Float32}(n)

@auto_hash_equals struct DBProperty
    name::String
    invert::Bool

    function DBProperty(n::DBNode)
        if startswith(n.key, "!")
            return new(n.key[2:end], true)
        else
            return new(n.key, false)
        end
    end
end

abstract type DBVariantAction end

@auto_hash_equals struct DBVariantDelete <: DBVariantAction
    atoms::Array{String}

    function DBVariantDelete(n::DBNode)
        new([a.key for a in n.value])
    end
end

@auto_hash_equals struct DBVariantRename <: DBVariantAction
    atoms::Dict{String, String}

    function DBVariantRename(n::DBNode)
        new(Dict(v.key => v.value for v in n.value))
    end
end

@auto_hash_equals struct DBVariant
    name::String

    actions::Array{DBVariantAction}
    properties::Array{DBProperty}

    function DBVariant(n::DBNode)
        name = n.key
        actions = []
        properties = []

        # variants can be of type delete or rename, and can carry properties
        for child in n.value
            if child.key == "Properties"
                properties = map(DBProperty, child.value)
            elseif child.key == "Delete"
                push!(actions, DBVariantDelete(child))
            elseif child.key == "Rename"
                push!(actions, DBVariantRename(child))
            else
                throw(ArgumentError("DBVariant: invalid format!"))
            end
        end

        new(name, actions, properties)
    end
end

@auto_hash_equals struct DBFragment{T<:Real}
    name::String
    path::String

    names::Array{String}
    atoms::Array{DBAtom{T}}
    bonds::Array{DBBond{T}}
    connections::Array{DBConnection{T}}
    properties::Array{DBProperty}
    variants::Array{DBVariant}
    
    function DBFragment{T}(n::DBNode) where {T<:Real}
        if startswith(n.key, "#include:")
            name = split(n.key, ":")[2]
            path = ball_data_path(split(n.value, ":")[1] * ".json")

            raw_fragment_data = JSON3.read(read(path, String), Array{DBNode})

            raw_names = filter(n -> n.key == "Names", raw_fragment_data)
            names = length(raw_names) == 1 ? [n.key for n in raw_names[1].value] : []
            
            raw_atoms = filter(n -> n.key == "Atoms", raw_fragment_data)
            atoms = length(raw_atoms) == 1 ? map(DBAtom{T}, raw_atoms[1].value) : []

            raw_bonds = filter(n -> n.key == "Bonds", raw_fragment_data)
            bonds = length(raw_bonds) == 1 ? map(b -> DBBond{T}(b), raw_bonds[1].value) : []

            raw_connections = filter(n -> n.key == "Connections", raw_fragment_data)
            connections = length(raw_connections) == 1 ? map(c -> DBConnection{T}(c), raw_connections[1].value) : []
            
            raw_properties = filter(n -> n.key == "Properties", raw_fragment_data)
            properties = length(raw_properties) == 1 ? map(DBProperty, raw_properties[1].value) : []

            raw_variants = filter(n -> n.key == "Variants", raw_fragment_data)
            variants = length(raw_variants) == 1 ? map(DBVariant, raw_variants[1].value) : []
            
            return new(name, path, names, atoms, bonds, connections, properties, variants)
        end

        throw(ArgumentError("DBFragment: invalid format!"))
    end
end
StructTypes.StructType(::Type{DBFragment}) = StructTypes.CustomStruct()

DBFragment(n::DBNode) = DBFragment{Float32}(n)

@auto_hash_equals struct DBNameMapping
    name::String
    maps_to::String

    mappings::Dict{String, String}

    function DBNameMapping(n::DBNode)
        if !startswith(n.key, "#include:")
            throw(ArgumentError("DBNameMapping: invalid format!"))
        end
    
        name = split(n.key, ":")[2]
        path = ball_data_path(split(n.value, ":")[1] * ".json")

        raw_mapping_data = JSON3.read(read(path, String), Array{DBNode})

        # the first node in the value list contains the reference
        if length(raw_mapping_data) == 0
            throw(ArgumentError("DBNameMapping: invalid format!"))
        end

        maps_to = raw_mapping_data[1].key

        mappings = Dict{String, String}(
            n.key => n.value for n in raw_mapping_data[2:end]
        )

        new(name, maps_to, mappings)
    end
end

@auto_hash_equals struct FragmentDB{T<:Real}
    fragments::OrderedDict{String, DBFragment{T}}
    name_mappings::OrderedDict{String, DBNameMapping}
    defaults::OrderedDict{String, String}

    function FragmentDB{T}(nodes::Vector{DBNode}) where {T<:Real}
        if length(nodes) == 3
            raw_fragments = filter(n -> n.key == "Fragments", nodes)
            raw_names     = filter(n -> n.key == "Names",     nodes)
            raw_defaults  = filter(n -> n.key == "Defaults",  nodes)

            if length(raw_fragments) == length(raw_names) == length(raw_defaults) == 1
                fragments = OrderedDict{String, DBFragment{T}}(
                    f.name => f for f in map(DBFragment{T}, raw_fragments[1].value)
                )
                
                name_mappings = OrderedDict{String, DBNameMapping}(
                    nm.name => nm for nm in map(DBNameMapping, raw_names[1].value)
                )

                defaults = OrderedDict{String, String}(
                    d.key => d.value for d in raw_defaults[1].value
                )

                return new(fragments, name_mappings, defaults)
            end
        end
        
        throw(ArgumentError("FragmentDB: invalid format!"))
    end

    function FragmentDB{T}(path::String = ball_data_path("fragments/Fragments.db.json")) where {T<:Real}
        jstring = read(path, String)
        
        JSON3.read(jstring, FragmentDB{T})
    end
end
StructTypes.StructType(::Type{FragmentDB{T}}) where {T} = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{FragmentDB{T}}) where {T} = Array{DBNode}

FragmentDB(path::String = ball_data_path("fragments/Fragments.db.json")) = FragmentDB{Float32}(path)


Base.show(io::IO, fdb::FragmentDB) = 
    print(io, 
        "FragmentDB with $(length(fdb.fragments)) fragments, " *
        "$(length(fdb.name_mappings)) mappings, $(length(fdb.defaults)) defaults.")