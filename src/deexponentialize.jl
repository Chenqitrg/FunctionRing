function Lie_action(ϵ::Vector{E}, f::HoloPoly{R,E}; trunc=O{R}(Inf)) where {R,E} # TODO: test
    f_temp = zero(R, E)
    for (n, ϵn) in enumerate(ϵ)
        f_temp = f_temp + ϵn * Witt_action(n, f)
    end
    return f_temp + trunc * f
end

function exp_Lie_action(ϵ::Vector{E}, f::HoloPoly{R,E}; trunc=O{R}(Inf)) where {R,E}
    f_temp_add = f
    f_temp_action = f
    n = 1
    while n < trunc
        f_temp_action = 1/n * Lie_action(ϵ, f_temp_action; trunc)
        f_temp_add = f_temp_add + f_temp_action
        n += 1
    end
    return f_temp_add
end

function deexponentialize(f::HoloPoly{Int,E}, trunc::O{Int}) where {E} # TODO: test
    if !(f[1] ≈ 1)
        throw(ArgumentError("Wrong polynomial"))
    end
    if trunc == O{Int}(Inf) || f.trunc == O{Int}(Inf)
        throw(ArgumentError("The truncation needs to be non-infinity"))
    end

    ϵ = [-f[2]]
    println("The coefficient of the 1-th generator has been calculated")
    poly = [HoloPoly{Int,E}([2], [f[2]])]
    real_trunc = min(floor(Int, trunc.trunc), floor(Int, f.trunc.trunc - 1))
    for n in 2:real_trunc
        poly_last = poly[end]
        for (i, p) in enumerate(poly)
            poly[i] = 1 / E(n - i + 1) * Lie_action(ϵ[1:i], p)
        end
        p_tot_coeff_np1 = sum(poly)[n+1]
        ϵn = p_tot_coeff_np1 - f[n+1]
        push!(ϵ, ϵn)
        push!(poly, poly_last + HoloPoly{Int,E}([n + 1], [-ϵn]))
        println("The coefficient of the $n-th generator has been calculated")
    end
    println("Final truncation: $real_trunc")
    return ϵ, O{Int}(real_trunc)
end