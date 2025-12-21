using ConcurrentCollections

struct VertexCache{N}
    storage::ConcurrentDict{CartesianIndex{N}, Tuple{Float64, SVector{N, Float64}}}
    tf::TestFunction
    lb::SVector{N, Float64}
    ub::SVector{N, Float64}
    cell_width::SVector{N, Float64}
end

function VertexCache(tf::TestFunction, divisions::NTuple{N, Int}) where N
    lb = SVector{N}(lb(tf))
    ub = SVector{N}(ub(tf))
    cell_width = (ub - lb) ./ SVector{N}(divisions)
    VertexCache(ConcurrentDict{CartesianIndex{N}, Tuple{Float64, SVector{N, Float64}}}(),
                tf, lb, ub, cell_width)
end

function get_vertex!(cache::VertexCache{N}, idx::CartesianIndex{N}) where N
    get!(cache.storage, idx) do
        x = cache.lb .+ (SVector(idx.I...) .- 1) .* cache.cell_width
        f = cache.tf.f(x)
        g = cache.tf.grad(x)
        (f, g)
    end
end