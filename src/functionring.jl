struct O{R}
    trunc::Float64
end

function Base.:+(f::O{R}, g::O{R}) where {R}
    if f.trunc == Inf && g.trunc == Inf
        return O{R}(Inf)
    elseif f.trunc == Inf
        return O{R}(g.trunc)
    elseif g.trunc == Inf
        return O{R}(f.trunc)
    else
        return O{R}(minimum([f.trunc, g.trunc]))
    end
end

function Base.:<(a, f::O{R}) where {R}
    return R(a) < f.trunc
end

function Base.:>(a, f::O{R}) where {R}
    return R(a) > f.trunc
end

function Base.:*(n::Real, f::O{R}) where {R}
    return O{R}(f.trunc + R(n))
end

function Base.:*(v::Vector, f::O{R}) where {R}
    if length(v) == 0
        return f
    else
        raised_pow = R(minimum(v)) + f.trunc
        return O{R}(raised_pow)
    end
end

struct HoloPoly{R, E}
    power_vect::Vector{R}
    coeff::Vector{E}
    trunc::O{R}

    function HoloPoly{R, E}(power_vect::Vector{R}, coeff::Vector{E}) where {R, E}
        new{R, E}(power_vect, coeff, O{R}(Inf))
    end
    function HoloPoly{R, E}(power_vect::Vector{R}, coeff::Vector{E}, trunc::Float64) where {R, E}
        new{R, E}(power_vect, coeff, O{R}(trunc))
    end
    function HoloPoly{R, E}(power_vect::Vector{R}, coeff::Vector{E}, trunc::O{R}) where {R, E}
        new{R, E}(power_vect, coeff, trunc)
    end
end

function Base.getindex(f::HoloPoly{R, E}, i::R) where {R, E}
    if i > f.trunc
        throw(ArgumentError("No higher element"))
    else
        index = findfirst(x -> x == i, f.power_vect)
        if index === nothing
            return 0
        else
            return f.coeff[index]
        end
    end
end

function Base.setindex!(f::HoloPoly{R, E}, setcoeff::E, n::R) where {R, E}
    if n in f.power_vect
        i = findfirst(x -> x == n, f.power_vect)
        f.coeff[i] = setcoeff
    elseif n < f.trunc
        push!(f.power_vect, n)
        push!(f.power_vect, setcoeff)
    end
end

function power_increasing(f::HoloPoly{R, E}) where {R, E}
    perm = sortperm(f.power_vect)
    new_pow = f.power_vect[perm]
    new_coeff = f.coeff[perm]
    trunc_pow = filter(x -> x < f.trunc, new_pow)
    return HoloPoly{R, E}(trunc_pow, new_coeff[1:length(trunc_pow)], f.trunc)
end

function zero_filter(f::HoloPoly{R, E}; tol = 1e-14) where {R, E}
    indices = findall(x -> abs(x) > tol, f.coeff)
    return HoloPoly{R, E}(f.power_vect[indices], f.coeff[indices], f.trunc)
end

function Base.:+(f::HoloPoly{R, E}, g::HoloPoly{R, E}) where {R, E}
    new_trunc = f.trunc + g.trunc
    new_power_vect = filter(x -> x < new_trunc, union(f.power_vect, g.power_vect))
    new_coeff = zeros(E, length(new_power_vect))
    for (i, n) in enumerate(new_power_vect)
       new_coeff[i] = f[n] + g[n]
    end
    return zero_filter(power_increasing(HoloPoly{R, E}(new_power_vect, new_coeff, new_trunc)))
end

function Base.:-(f::HoloPoly{R, E}, g::HoloPoly{R, E}) where {R, E}
    new_power_vect = union(f.power_vect, g.power_vect)
    new_coeff = zeros(E, length(new_power_vect))
    for (i, n) in enumerate(new_power_vect)
       new_coeff[i] = f[n] - g[n]
    end
    return zero_filter(power_increasing(HoloPoly{R, E}(new_power_vect, new_coeff, f.trunc + g.trunc)))
end

function Base.:-(f::HoloPoly{R, E}) where {R, E}
    return zero_filter(power_increasing(HoloPoly{R, E}(f.power_vect, -f.coeff, f.trunc)))
end

function Base.:*(f::HoloPoly{R, E}, g::HoloPoly{R, E}) where {R, E}
    fg_power_vect = R[]
    fg_trunc = f.power_vect * g.trunc + g.power_vect * f.trunc

    for n_f in f.power_vect, n_g in g.power_vect
        temp_n = n_f + n_g
        if !(temp_n in fg_power_vect) && temp_n < fg_trunc
            push!(fg_power_vect, temp_n)
        end
    end
    fg_coeff = zeros(E, length(fg_power_vect))
    for (i_f, n_f) in enumerate(f.power_vect), (i_g, n_g) in enumerate(g.power_vect)
        n_fg = n_f + n_g
        if n_fg < fg_trunc
            i_fg = findfirst(x -> x == n_fg, fg_power_vect)
            fg_coeff[i_fg] += f.coeff[i_f] * g.coeff[i_g]
        else
            continue
        end
    end
    return zero_filter(power_increasing(HoloPoly{R, E}(fg_power_vect, fg_coeff, fg_trunc)))
end

function Base.:*(a::E, f::HoloPoly{R, E}) where {R, E}
    return zero_filter(power_increasing(HoloPoly{R, E}(f.power_vect, a * f.coeff, f.trunc)))
end

function Base.:^(f::HoloPoly{R, E}, n::Int) where {R, E}
    if n == 0
        return HoloPoly{R, E}([0], [E(1)])
    else
        running_f = f
        power = 1
        while power < n
            running_f *= f
            power += 1
        end
        return running_f
    end
end

function evaluation(f::HoloPoly{R, E}, z::E) where {R, E}
    running_numb = E(0)
    for (i, n) in enumerate(f.power_vect)
        running_numb += z^n * f.coeff[i]
    end
    return running_numb
end

f = HoloPoly{Int, Float64}([0, 2, 1, 3, 4], [1.0, 2.0, 1.2, 0.5,3.0])
# # power_increasing(f)

g = HoloPoly{Int, Float64}([0, 1, 3, 5], [1.5, 2.0, 4.2, 3.0], O{Int}(8.0))

# # 0.0 * f

# f + g

g * g * g

g^3

evaluation(f, 0.5)


# (f - g) * g * f
# O{Int}(3.0)

# HoloPoly{Int64, Float64}([0, 2, 1], ComplexF64[2.5 + 0.0im, 2.0 + 0.0im, 3.2 + 0.0im], O{Int64}(3.0))
# function Base.:*(a::R, f::HoloPoly{R}) where {R}
#     newf = copy(f)

# end