using StaticArrays
using Base.Threads

struct VertexCache{N}
    storage::Dict{CartesianIndex{N}, Tuple{Float64, SVector{N, Float64}}}
    lock::ReentrantLock
    tf::TestFunction
    lb::SVector{N, Float64}
    ub::SVector{N, Float64}
    cell_width::SVector{N, Float64}
end

function VertexCache(tf::TestFunction, divisions::NTuple{N, Int}) where N
    lb = SVector{N}(lb(tf))
    ub = SVector{N}(ub(tf))
    cell_width = (ub - lb) ./ SVector{N}(divisions)
    VertexCache(
        Dict(),
        ReentrantLock(),
        tf, lb, ub, cell_width
    )
end

function get_vertex!(cache::VertexCache{N}, idx::CartesianIndex{N}) where N
    lock(cache.lock) do
        get!(cache.storage, idx) do
            x = cache.lb .+ (SVector(idx.I...) .- 1) .* cache.cell_width
            f = cache.tf.f(x)
            g = cache.tf.grad(x)
            (f, g)
        end
    end
end