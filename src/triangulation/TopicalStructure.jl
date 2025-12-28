module TopicalStructure
using StaticArrays, DataStructures

export TopicalManager, Vertex, Simplex, add_vertex!, add_simplex!, all_simplex_ids

struct Vertex{N}
    id::Int
    pos::SVector{N, Float64}
    val::Float64
end

struct Simplex{N}
    id::Int
    vertices::Vector{Int}
    min_val::Float64
end

mutable struct TopicalManager{N}
    vertices::Dict{Int, Vertex{N}}
    simplices::Dict{Int, Simplex{N}}
    star_map::Dict{Int, Vector{Int}}
    
    function TopicalManager{N}() where N
        new{N}(Dict{Int, Vertex{N}}(), Dict{Int, Simplex{N}}(), Dict{Int, Vector{Int}}())
    end
end

function add_vertex!(tm::TopicalManager{N}, pos::SVector{N, Float64}, val::Float64) where N
    v_id = length(tm.vertices) + 1
    v = Vertex{N}(v_id, pos, val)
    tm.vertices[v_id] = v
    tm.star_map[v_id] = Int[]
    return v
end

function add_simplex!(tm::TopicalManager{N}, v_ids::Vector{Int}) where N
    s_id = length(tm.simplices) + 1
    min_v = minimum(tm.vertices[vid].val for vid in v_ids)
    s = Simplex{N}(s_id, v_ids, min_v)
    tm.simplices[s_id] = s
    for vid in v_ids
        push!(tm.star_map[vid], s_id)
    end
    return s_id
end

function all_simplex_ids(tm::TopicalManager)
    return collect(keys(tm.simplices))
end

end # module