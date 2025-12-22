function Lie_action(ϵ::Vector{E}, f::HoloPoly{R,E}) where {R,E} # TODO: test
    f_temp = zero(R, E)
    for (n, ϵn) in enumerate(ϵ)
        f_temp = f_temp + ϵn * Witt_action(n, f)
    end
    return f_temp
end

function deexponentialize(f::HoloPoly{Int,E}, trunc::Int) where {E} # TODO: test
    ϵ = [-f[2]]
    poly = [HoloPoly{Int,E}([2], [f[2]])]
    for n in 2:trunc
        @show n, ϵ
        poly_last = poly[end]
        for (i, p) in enumerate(poly)
            poly[i] = 1 / E(n - i + 1) * Lie_action(ϵ[1:i], p)
        end
        p_tot_coeff_np1 = sum(poly)[n+1]
        ϵn = p_tot_coeff_np1 - f[n+1]
        push!(ϵ, ϵn)
        push!(poly, poly_last + HoloPoly{Int,E}([n + 1], [-ϵn]))
    end
    return ϵ, poly
end