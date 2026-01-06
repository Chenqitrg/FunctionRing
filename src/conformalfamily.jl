level(v::Vector{Int}) = -sum(v)
left_action(n::Int, v::Vector{Int}) = [n; v...]
iscanonical(v::Vector{Int}) = isempty(v) || (issorted(v; lt=(<)) && v[end] < 0)
canonicalbasis(n::Integer) = n == 0 ? (Int[],) : (.-p for p in partitions(n))
dim_level(level::Integer) = length(partitions(level))

function iscanonical(bs::Vector{Vector{Int}})
    level_first = level(bs[1])
    return all(level(x) == level_first && iscanonical(x) for x in bs)
end

function to_conformal_index(b::Vector{Int})
    level_b = level(b)
    i = 1
    for p in canonicalbasis(level_b)
        p == b && return level_b, i
        i += 1
    end
    return nothing
end

struct Descendent{E}
    basis::Vector{Vector{Int}}
    vector::Vector{E}
end

function Descendent{E}(level::Int, vector::SparseVector{E,Int64}) where {E}
    I = vector.nzind
    vals = vector.nzval
    basis = Vector{Int}[]
    for (i, p) in enumerate(canonicalbasis(level))
        (i in I) && push!(basis, p)
    end
    return Descendent{E}(basis[I], vals)
end

function Base.:(==)(f::Descendent{E}, g::Descendent{E}) where {E}
    return f.basis == g.basis && f.vector == g.vector
end

function Base.getindex(f::Descendent{E}, v::Vector{Int}) where {E}
    i = findfirst(==(v), f.basis)
    return i === nothing ? zero(E) : f.vector[i]
end

function Base.:+(f::Descendent{E}, g::Descendent{E}) where {E}
    all_basis = union(f.basis, g.basis)
    new_vector = [f[v] + g[v] for v in all_basis]
    return Descendent{E}(all_basis, new_vector)
end

function Base.:*(a, g::Descendent{E}) where {E}
    return Descendent{E}(g.basis, E(a) * g.vector)
end

function iscanonical(f::Descendent{E}) where {E}
    if length(f.basis) != length(f.vector)
        return false
    elseif length(f.basis) == 0
        return true
    elseif length(f.basis) != length(unique(f.basis))
        return false
    end
    level_f = level(f.basis[1])
    return all(iscanonical(b) && level(b) == level_f for b in f.basis)
end

function once_canonicalize(b::Vector{Int})
    !(b[1] < 0 && iscanonical(b[2:end])) && throw(ArgumentError("The vector $b is not suitable for once_canonicalize"))
    iscanonical(b) && return Descendent([b], [1.0])
    first_level = b[1]
    second_level = b[2]
    newbasis1 = [first_level + second_level; b[3:end]...]
    vector = Descendent([newbasis1], [Float64(first_level - second_level)])

    newbasis2 = [first_level; b[3:end]...]
    shuffle1 = once_canonicalize(newbasis2)
    for basis2 in shuffle1.basis
        shuffle2 = shuffle1[basis2] * once_canonicalize(left_action(second_level, basis2))
        vector = vector + shuffle2
    end
    return vector
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, row::Vector{Int}, col::Vector{Int}) where {E}
    rowind, colind = to_conformal_index(row), to_conformal_index(col)
    return mat[rowind[2], colind[2]]
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, row::Vector{Int}, ::Colon) where {E}
    rowind = to_conformal_index(row)
    return mat[rowind[2], :]
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, ::Colon, col::Vector{Int}) where {E}
    colind = to_conformal_index(col)
    return mat[:, colind[2]]
end

function Base.getindex(vec::SparseVector{E,Int64}, ind::Vector{Int}) where {E}
    i = to_conformal_index(ind)
    return vec[i[2]]
end

function Base.setindex!(mat::SparseMatrixCSC{E,Int64}, ele::E, row::Vector{Int}, col::Vector{Int}) where {E}
    rowind = to_conformal_index(row)
    colind = to_conformal_index(col)
    mat[rowind[2], colind[2]] = ele
end

function Base.setindex!(vec::SparseVector{E,Int64}, ele::E, ind::Vector{Int}) where {E}
    i = to_conformal_index(ind)
    vec[i[2]] = ele
end

function vira_generator_level_col_solver(h::Float64, c::Float64, col::Vector{Int}, operator_level::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    col_level = level(col)
    (col_level == 1 && operator_level == 1) && return sparse([2 * h])
    row_level = col_level - operator_level
    row_dim = dim_level(row_level)

    highest_operator = col[1]
    remaining_col = col[2:end]
    remaining_level = col_level + highest_operator
    new_leading_operator = operator_level + highest_operator

    if new_leading_operator > 0
        vector = (operator_level - highest_operator) * reps[remaining_level][new_leading_operator][:, remaining_col]
    elseif new_leading_operator == 0
        vector = sparsevec([to_conformal_index(remaining_col)[2]], (operator_level - highest_operator) * (h + remaining_level) + c / 12 * (operator_level^3 - operator_level), row_dim)
    else
        canonicalized_basis = once_canonicalize(left_action(new_leading_operator, remaining_col))
        inds = Int[]
        vals = Float64[]
        for row_vec in canonicalized_basis.basis
            row_val = (operator_level - highest_operator) * canonicalized_basis[row_vec]
            if row_val != 0
                push!(inds, to_conformal_index(row_vec)[2])
                push!(vals, row_val)
            end
        end
        vector = sparsevec(inds, vals, row_dim)
    end

    remaining_shifted_level = remaining_level - operator_level
    remaining_shifted_level < 0 && return vector

    remaining_action_mat = reps[remaining_level][operator_level]
    j = to_conformal_index(remaining_col)[2]
    i_s = remaining_action_mat.rowval[remaining_action_mat.colptr[j] : remaining_action_mat.colptr[j+1] - 1]
    for (i,reduced_row) in enumerate(canonicalbasis(remaining_shifted_level))
        !(i in i_s) && continue
        second_canonicalized_basis = remaining_action_mat[i, j] * once_canonicalize(left_action(highest_operator, reduced_row))
        for row in second_canonicalized_basis.basis
            vector[row] += second_canonicalized_basis[row]
        end
    end
    return vector
end

function vira_generator_level_solver(h::Float64, c::Float64, col_level::Int, operator_level::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    colptr = Int[1]
    rowind = Int[]
    nzval = Float64[]

    nnz = 0
    for col in canonicalbasis(col_level)
        v = vira_generator_level_col_solver(h, c, col, operator_level, reps)
        # println("$col has been obtained")
        append!(rowind, v.nzind)
        append!(nzval, v.nzval)

        nnz += length(v.nzind)
        push!(colptr, nnz + 1)
    end

    nrow = length(partitions(col_level - operator_level))
    ncol = length(colptr) - 1

    M = SparseMatrixCSC(nrow, ncol, colptr, rowind, nzval)
    return M
end

function vira_level_solver(h::Float64, c::Float64, level_number::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    matrix_vector = SparseMatrixCSC{Float64,Int64}[]
    for vira_n in 1:level_number
        @time matrix = vira_generator_level_solver(h, c, level_number, vira_n, reps)
        println("level $level_number operator $vira_n has been calculated")
        push!(matrix_vector, matrix)
    end
    return matrix_vector
end

function vira_iter_solver(h::Float64, c::Float64, trunc_level::Int)
    reps = Vector{SparseMatrixCSC{Float64,Int64}}[]
    for level_number in 1:trunc_level
        matrices = vira_level_solver(h, c, level_number, reps)
        push!(reps, matrices)
    end
    return reps
end

function _vira_block_embedding(rep_mats::Vector{SparseMatrixCSC{Float64,Int64}}, operator_level::Int, trunc_level::Int)
    dim_row = sum(dim_level(level) for level in 0:trunc_level)
    dim_col = dim_row
    block_mat = spzeros(Float64, dim_row, dim_col)
    col_offset = sum(dim_level(level) for level in 0:operator_level-1)
    row_offset = 0
    for col_level in operator_level:trunc_level
        row_level = col_level - operator_level
        row_dim = dim_level(row_level)
        col_dim = dim_level(col_level)
        block_mat[row_offset.+(1:row_dim), col_offset.+(1:col_dim)] = rep_mats[col_level-operator_level+1]
        col_offset += col_dim
        row_offset += row_dim
    end
    return block_mat
end

function vira_reps_total(reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    trunc_level = length(reps)
    block_mats = SparseMatrixCSC{Float64,Int64}[]
    for level in 1:trunc_level
        mats = [reps[n][level] for n in level:trunc_level]
        block_mat = _vira_block_embedding(mats, level, trunc_level)
        push!(block_mats, block_mat)
    end
    return block_mats
end

function gram_matrix_level(level::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}}, gram_lower::Vector{SparseMatrixCSC{Float64,Int64}})
    dim = dim_level(level)
    gram = spzeros(Float64, dim, dim)
    for bra in canonicalbasis(level), ket in canonicalbasis(level)
        level == 0 && (gram[bra, ket] = 1.0; continue)
        highest_level_bra = bra[1]
        remaining_bra = bra[2:end]
        remaining_bra_level = level + highest_level_bra
        actioning_operator = -highest_level_bra

        gram[bra, ket] = dot(gram_lower[remaining_bra_level+1][remaining_bra, :], reps[level][actioning_operator][:, ket])
    end
    return gram
end

function gram_matrix(reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    gram_lower = SparseMatrixCSC{Float64,Int64}[]
    for level in 0:length(reps)
        gram = gram_matrix_level(level, reps, gram_lower)
        push!(gram_lower, gram)
    end
    return gram_lower
end