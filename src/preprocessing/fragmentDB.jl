using StructTypes
using JSON3
using AutoHashEquals

import DataStructures: OrderedCollections.OrderedDict

export FragmentDB, get_reference_fragment

@auto_hash_equals struct DBNode
    key::String
    value::Union{AbstractArray{DBNode}, String}
end
StructTypes.StructType(::Type{DBNode}) = StructTypes.Struct()

@auto_hash_equals struct DBAtom{T<:Real}
    name::String
    element::ElementType
    r::Vector3

    function DBAtom{T}(name::String, element::ElementType, r::Vector3) where {T<:Real}
        new(name, element, r)
    end

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

@auto_hash_equals struct DBBond
    number::Int
    a1::String
    a2::String
    order::BondOrderType

    function DBBond(number::Int, a1::String, a2::String, order::BondOrderType)
        new(number, a1, a2, order)
    end

    function DBBond(n::DBNode)
        number = parse(Int, n.key)

        raw_data = split(n.value)

        a1 = raw_data[1]
        a2 = raw_data[2]
        
        order = to_bond_order(raw_data[3])

        new(number, a1, a2, order)
    end
end

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
    value::Bool

    function DBProperty(n::DBNode)
        if startswith(n.key, "!")
            return new(n.key[2:end], false)
        else
            return new(n.key, true)
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

@auto_hash_equals struct DBVariant{T<:Real}
    name::String

    atoms::Array{DBAtom{T}}
    bonds::Array{DBBond}

    actions::Array{DBVariantAction}
    properties::Array{DBProperty}

    function DBVariant{T}(n::DBNode, atoms::Array{DBAtom{T}}, bonds::Array{DBBond}) where {T<:Real}
        name = n.key
        actions = []
        properties = []

        # variants can be of type delete or rename, and can carry properties
        for child in n.value
            if child.key == "Properties"
                properties = map(DBProperty, child.value)
            elseif child.key == "Delete"
                db_delete = DBVariantDelete(child)

                atoms = filter(a -> a.name ∉ db_delete.atoms, atoms)
                bonds = filter(b -> b.a1 ∉ db_delete.atoms && b.a2 ∉ db_delete.atoms, bonds)
                
                push!(actions, db_delete)
            elseif child.key == "Rename"
                db_rename = DBVariantRename(child)

                atoms = map(
                    a -> DBAtom{T}(get(db_rename.atoms, a.name, a.name), a.element, a.r), atoms)
                bonds = map(
                    b -> DBBond(b.number,
                            get(db_rename.atoms, b.a1, b.a1),
                            get(db_rename.atoms, b.a2, b.a2),
                            b.order), bonds)
                                
                push!(actions, db_rename)
            else
                throw(ArgumentError("DBVariant: invalid format!"))
            end
        end

        new(name, atoms, bonds, actions, properties)
    end
end

@auto_hash_equals struct DBFragment{T<:Real}
    name::String
    path::String

    names::Array{String}
    atoms::Array{DBAtom{T}}
    bonds::Array{DBBond}
    connections::Array{DBConnection{T}}
    properties::Array{DBProperty}
    variants::Array{DBVariant{T}}
    
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
            bonds = length(raw_bonds) == 1 ? map(b -> DBBond(b), raw_bonds[1].value) : []

            raw_connections = filter(n -> n.key == "Connections", raw_fragment_data)
            connections = length(raw_connections) == 1 ? map(c -> DBConnection{T}(c), raw_connections[1].value) : []
            
            raw_properties = filter(n -> n.key == "Properties", raw_fragment_data)
            properties = length(raw_properties) == 1 ? map(DBProperty, raw_properties[1].value) : []

            raw_variants = filter(n -> n.key == "Variants", raw_fragment_data)
            variants = length(raw_variants) == 1 ? map(v-> DBVariant{T}(v, atoms, bonds), raw_variants[1].value) : []
            
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

@auto_hash_equals struct FragmentDB
    fragments::OrderedDict{String, DBFragment}
    name_mappings::OrderedDict{String, DBNameMapping}
    defaults::OrderedDict{String, String}

    function FragmentDB(nodes::Vector{DBNode})
        if length(nodes) == 3
            raw_fragments = filter(n -> n.key == "Fragments", nodes)
            raw_names     = filter(n -> n.key == "Names",     nodes)
            raw_defaults  = filter(n -> n.key == "Defaults",  nodes)

            if length(raw_fragments) == length(raw_names) == length(raw_defaults) == 1
                fragments = OrderedDict{String, DBFragment}(
                    f.name => f for f in map(DBFragment, raw_fragments[1].value)
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

    function FragmentDB(path::String = ball_data_path("fragments/fragments.db.json"))
        jstring = read(path, String)
        
        JSON3.read(jstring, FragmentDB)
    end
end
StructTypes.StructType(::Type{FragmentDB}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{FragmentDB}) = Array{DBNode}

function get_reference_fragment(f::Fragment{T}, fdb::FragmentDB) where {T<:Real}
    # first, try to find the fragment in the database
    if f.name ∉ keys(fdb.fragments)
        return nothing
    end

    db_fragment = fdb.fragments[f.name]

    # does the fragment have variants?
    if length(db_fragment.variants) == 1
        return db_fragment.variants[1]
    end

    # now, find the variant that best matches the fragment
    # This returns N/C terminal variants for fragments that
    # have the corresponding properties set or cystein variants
    # without thiol hydrogen if the disulphide bond property
    # is set

    # First, check for two special properties of amino acids:
    # C_TERMINAL and N_TERMINAL
    # They are usually not set, so set them here
    if is_c_terminal(f)
        f.properties["C_TERMINAL"] = true
    end

    if is_n_terminal(f)
        f.properties["N_TERMINAL"] = true
    end

    if is_3_prime(f)
        f.properties["3_PRIME"] = true
    end

    if is_5_prime(f)
        f.properties["5_PRIME"] = true
    end

    # the number of properties that matched
    # the fragment with the largest number of matched
    # properties is returned
    number_of_properties = -1
    property_difference = -1
    best_number_of_properties = -1
    best_property_difference = 10000

    best_variant = nothing

    # iterate over all variants of the fragment and compare the properties
    for var in db_fragment.variants
        var_props = Dict(p.name => p.value for p in var.properties)

        # count the difference in the number of set properties
        property_difference = abs(length(findall(identity, f.properties)) - length(findall(identity, var_props)))

        # and count the properties fragment and variant have in common
        number_of_properties = length(f.properties ∩ var_props)

        @debug "Considering variant $(var.name). # properties: $(number_of_properties)"

        if ((number_of_properties > best_number_of_properties)
            || (   (number_of_properties == best_number_of_properties)
                && (property_difference < best_property_difference)))
            best_variant = var
            best_number_of_properties = number_of_properties
            best_property_difference  = property_difference
        end
    end

    best_variant
end

Base.show(io::IO, fdb::FragmentDB) = 
    print(io, 
        "FragmentDB with $(length(fdb.fragments)) fragments, " *
        "$(length(fdb.name_mappings)) mappings, $(length(fdb.defaults)) defaults.")