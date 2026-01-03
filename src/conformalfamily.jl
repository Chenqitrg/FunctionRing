struct DescendentBasis
    basis::Vector{Pair{Int,Int}}
    function DescendentBasis(v::Vector{Pair{Int,Int}})
        return new(filter(x -> x[2] != 0, v))
    end
end

function level(v::Vector{Pair{Int,Int}})
    return isempty(v) ? 0 : sum(-i * n for (i, n) in v)
end
level(v::DescendentBasis) = level(v.basis)

Base.:(==)(u::DescendentBasis, v::DescendentBasis) = u.basis == v.basis

function left_action(n::Int, b::DescendentBasis)
    if isempty(b.basis)
        return DescendentBasis([n => 1])
    end
    leading_level, leading_power = b.basis[1]
    if leading_level == n
        return DescendentBasis([leading_level => leading_power + 1; b.basis[2:end]...])
    else
        return DescendentBasis([n => 1; b.basis...])
    end
end

function iscanonical(b::DescendentBasis)
    return isempty(b.basis) || (issorted(map(x -> x[1], b.basis); lt=(<)) && b.basis[end][1] < 0) && unique(map(x -> x[1], b.basis)) == map(x -> x[1], b.basis)
end

function iscanonical(bs::Vector{DescendentBasis})
    level_first = level(bs[1])
    return all(level(x) == level_first && iscanonical(x) for x in bs)
end

function _canonical_upperbound(n::Int, bound::Int)
    if bound == 1
        return [[-1 => n]]
    elseif n == 0 && bound == 0
        return Vector{Pair{Int,Int}}[]
    end
    basis = Vector{Pair{Int,Int}}[]
    for power in 0:n÷bound
        temp_pair = [-bound => power]
        remaining = n - power * bound
        if remaining == 0
            push!(basis, temp_pair)
        else
            remaining_results = _canonical_upperbound(remaining, minimum([bound - 1, remaining]))
            append!(basis, map(x -> [temp_pair; x...], remaining_results))
        end
    end
    return basis
end

function _to_vector(v::Vector{Pair{Int,Int}})
    level_v = level(v)
    vector_v = zeros(Int, level_v)
    if level_v == 0
        return vector_v
    end
    for n in -level_v:-1
        pair = findfirst(x -> x[1] == n, v)
        vector_v[level_v+n+1] = pair === nothing ? 0 : v[pair][2]
    end
    return vector_v
end

weight_key(v::Vector{Pair{Int,Int}}) = (level(v), _to_vector(v))

function canonicalbasis(n::Int)
    if n == 0
        return [DescendentBasis(Pair{Int,Int}[])]
    else
        return map(DescendentBasis, sort(_canonical_upperbound(n, n); by=weight_key, rev=true))
    end
end

canonicalbasis_plain(n::Int) = map(x -> x.basis, canonicalbasis(n))

function to_conformal_index(b::DescendentBasis)
    level_b = level(b)
    i = findfirst(==(b), CANONICALBASIS[level_b+1])
    return level_b, i
end

struct Descendent{E}
    basis::Vector{DescendentBasis}
    vector::Vector{E}
end

function Descendent{E}(basis::Vector{Vector{Pair{Int,Int}}}, vector) where {E}
    return Descendent{E}(map(DescendentBasis, basis), vector)
end

function Base.:(==)(f::Descendent{E}, g::Descendent{E}) where {E}
    return f.basis == g.basis && f.vector == g.vector
end

function Base.getindex(f::Descendent{E}, v::DescendentBasis) where {E}
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

function once_canonicalize(b::DescendentBasis)
    if iscanonical(b)
        return Descendent([b], [1.0])
    elseif b.basis[1][1] == b.basis[2][1]
        newbasis = DescendentBasis([b.basis[1][1] => b.basis[1][2] + b.basis[2][2]; b.basis[3:end]...])
        return Descendent([newbasis], [1.0])
    end
    first_level, first_power = b.basis[1]
    second_level, second_power = b.basis[2]
    newbasis1 = DescendentBasis([[first_level + second_level => 1, second_level => second_power - 1]; b.basis[3:end]...])
    vector = Descendent([newbasis1], [Float64(first_level - second_level)])

    newbasis2 = DescendentBasis([[first_level => 1, second_level => second_power - 1]; b.basis[3:end]...])
    shuffle1 = once_canonicalize(newbasis2)
    for basis2 in shuffle1.basis
        shuffle2 = shuffle1[basis2] * once_canonicalize(left_action(second_level, basis2))
        vector += shuffle2
    end
    return vector
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, row::DescendentBasis, col::DescendentBasis) where {E}
    rowind, colind = to_conformal_index(row), to_conformal_index(col)
    return mat[rowind[2], colind[2]]
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, row::DescendentBasis, ::Colon) where {E}
    rowind = to_conformal_index(row)
    return mat[rowind[2], :]
end

function Base.getindex(mat::SparseMatrixCSC{E,Int64}, ::Colon, col::DescendentBasis) where {E}
    colind = to_conformal_index(col)
    return mat[:, colind[2]]
end

function Base.setindex!(mat::SparseMatrixCSC{E,Int64}, ele::E, row::DescendentBasis, col::DescendentBasis) where {E}
    rowind = to_conformal_index(row)
    colind = to_conformal_index(col)
    mat[rowind[2], colind[2]] = ele
end

function vira_generator_level_solver(h::Float64, c::Float64, col_level::Int, operator_level::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    if col_level == 1 && operator_level == 1
        return sparse([2 * h;;])
    end
    row_level = col_level - operator_level
    row_basis = CANONICALBASIS[row_level+1]
    col_basis = CANONICALBASIS[col_level+1]
    matrix = spzeros(Float64, length(row_basis), length(col_basis))
    for col in col_basis
        highest_operator, highestpower = col.basis[1]
        remaining_col = DescendentBasis([highest_operator => highestpower - 1; col.basis[2:end]...])
        remaining_level = col_level + highest_operator
        new_leading_operator = operator_level + highest_operator

        if new_leading_operator > 0
            for row in row_basis
                matrix[row, col] = (operator_level - highest_operator) * reps[remaining_level][new_leading_operator][row, remaining_col]
            end
        elseif new_leading_operator == 0
            matrix[remaining_col, col] = (operator_level - highest_operator) * (h + remaining_level) + c / 12 * (operator_level^3 - operator_level)
        else
            canonicalized_basis = once_canonicalize(left_action(new_leading_operator, remaining_col))
            for row in row_basis
                matrix[row, col] = (operator_level - highest_operator) * canonicalized_basis[row]
            end
        end

        remaining_shifted_level = remaining_level - operator_level
        remaining_shifted_level < 0 && continue

        for reduced_row in CANONICALBASIS[remaining_shifted_level+1]
            if isempty(reduced_row.basis) || (reduced_row.basis[1][1] >= highest_operator)
                matrix[left_action(highest_operator, reduced_row), col] += reps[remaining_level][operator_level][reduced_row, remaining_col]
            end
        end
    end
    return matrix
end

function vira_level_solver(h::Float64, c::Float64, level_number::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}})
    matrix_vector = SparseMatrixCSC{Float64,Int64}[]
    for vira_n in 1:level_number
        matrix = vira_generator_level_solver(h, c, level_number, vira_n, reps)
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

function gram_matrix_level(level::Int, reps::Vector{Vector{SparseMatrixCSC{Float64,Int64}}}, gram_lower::Vector{SparseMatrixCSC{Float64,Int64}})
    basis = CANONICALBASIS[level+1]
    dim = length(basis)
    gram = spzeros(Float64, dim, dim)
    for bra in basis, ket in basis
        level == 0 && (gram[bra, ket] = 1.0; continue)
        highest_level_bra, highest_power_bra = bra.basis[1]
        remaining_bra = DescendentBasis([highest_level_bra => highest_power_bra - 1; bra.basis[2:end]...])
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