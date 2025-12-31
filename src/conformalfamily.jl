level(v::Vector{Int}) = isempty(v) ? 0 : sum(i * n for (i, n) in enumerate(v))
weight_key(v::Vector{Int}) = (level(v), v)

function partitions(n::Int)
    return _partitions(n, 1)
end

function _partitions(n::Int, minpart::Int)
    if n == 0
        return [Int[]]  # one partition: the empty partition
    end
    parts = Vector{Vector{Int}}()
    for i in minpart:n
        if i > n
            break
        end
        for tail in _partitions(n - i, i)
            push!(parts, [i; tail])
        end
    end
    return parts
end

function partition_to_basis(v::Vector{Int})
    basis_vec = Vector{Int}()
    if isempty(v)
        return basis_vec
    end
    for i in 1:v[end]
        push!(basis_vec, count(==(i), v))
    end
    return basis_vec
end

function canonicalbasis(n::Int)
    return map(partition_to_basis, partitions(n))
end

@views function tailview(v::Vector{Int})
    lastidx = findlast(!iszero, v)
    lastidx === nothing && return view(v, 1:0)
    return view(v, 1:lastidx)
end

struct Descendent{E}
    basis::Vector{Vector{Int}}
    vector::Vector{E}
    function Descendent{E}(basis::Vector{Vector{Int}}, vector::Vector{E}) where {E}
        perm = sortperm(basis; by=weight_key)
        return new(map(tailview, basis[perm]), vector[perm])
    end
end

function Base.:(==)(f::Descendent{E}, g::Descendent{E}) where {E}
    return f.basis == g.basis && f.vector == g.vector
end

function Base.:(≈)(f::Descendent{E}, g::Descendent{E}) where {E}
    return f.basis == g.basis && f.vector ≈ g.vector
end

function Base.getindex(f::Descendent{E}, v::Vector{Int}) where {E}
    i = findfirst(x -> x == tailview(v), f.basis)
    return i === nothing ? zero(E) : f.vector[i]
end

function Base.:+(f::Descendent{E}, g::Descendent{E}) where {E}
    all_basis = union(f.basis, g.basis)
    new_vector = [f[v] + g[v] for v in all_basis]
    return Descendent{E}(all_basis, new_vector)
end

function Base.:*(a::E, g::Descendent{E}) where {E}
    return Descendent{E}(g.basis, a * g.vector)
end

struct ConformalOperator{E}
    rowbasis::Vector{Vector{Int}}
    colbasis::Vector{Vector{Int}}
    matrix::Array{E,2}
    function ConformalOperator(row::Vector{Vector{Int}}, col::Vector{Vector{Int}}, matrix_ele::Array{E,2}) where {E}
        rowperm = sortperm(row; by=weight_key)
        colperm = sortperm(col; by=weight_key)
        new_rowbasis = map(tailview, row[rowperm])
        new_colbasis = map(tailview, col[colperm])
        new_matrix = matrix_ele[rowperm, colperm]
        return new{E}(new_rowbasis, new_colbasis, new_matrix)
    end
end