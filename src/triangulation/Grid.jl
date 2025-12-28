module Grid
using StaticArrays

export AbstractGrid, GridStructure, GridPoint, VertexCache

abstract type AbstractGrid end

struct GridPoint{N}
    idx::CartesianIndex{N}
    pos::SVector{N, Float64}
    v_id::Int
end

struct VertexCache{N}
    points::Dict{CartesianIndex{N}, GridPoint{N}}
    VertexCache{N}() where N = new{N}(Dict{CartesianIndex{N}, GridPoint{N}}())
end

struct GridStructure{N} <: AbstractGrid
    lower::SVector{N, Float64}
    upper::SVector{N, Float64}
    dims::SVector{N, Int}
    steps::SVector{N, Float64}
    cache::VertexCache{N}

    function GridStructure(lower::Vector{Float64}, upper::Vector{Float64}, dims::Vector{Int})
        N = length(lower)
        steps = SVector{N, Float64}([(upper[i] - lower[i]) / (dims[i] - 1) for i in 1:N]...)
        new{N}(SVector{N}(lower...), SVector{N}(upper...), SVector{N}(dims...), 
               steps, VertexCache{N}())
    end
end

function calculate_pos(grid::GridStructure{N}, idx::CartesianIndex{N}) where N
    return grid.lower .+ (SVector{N, Float64}(idx.I...) .- 1.0) .* grid.steps
end

end # module