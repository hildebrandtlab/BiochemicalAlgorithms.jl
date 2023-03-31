export AtomTypeTemplate

@auto_hash_equals struct AtomTypeTemplate{T<:Real}
    type_name::String
    charge::T
end