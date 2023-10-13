export IndexedDataFrame

struct IndexedDataFrame
    df::DataFrame
    index::Dict{Int, Int}

    function IndexedDataFrame(values::Vector{T}) where T
        df = DataFrame(values)
        index = Dict{Int, Int}()

        new(df, index)
    end
end

function Base.push!(idf::IndexedDataFrame, v::Tuple)
    new_rownumber = last(push!(idf.df, v))
    idf.index[v[1]] = rownumber(new_rownumber)
end

function Base.push!(idf::IndexedDataFrame, v::NamedTuple)
    new_rownumber = last(push!(idf.df, v))
    idf.index[v.idx] = rownumber(new_rownumber)
end

@inline function _row_by_idx(idf::IndexedDataFrame, idx::Int)
    idf.index[idx]
end