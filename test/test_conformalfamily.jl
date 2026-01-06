println("----------------")
println("Conformal family")
println("----------------")

@testset "To index" begin
    set = [[-2], [-1, -1]]
    b = [-2]
    @test findfirst(==(b), set) == 1
end

@testset "Union basis" begin
    b1 = [-3, -2]
    b2 = [-5]
    u1 = [b1]
    u2 = [b2, b1]
    @test union(u1, u2) == [[-3, -2], [-5]]
    v1 = [[-3, -2]]
    v2 = [[-5], [-3, -2]]
    @test union(v1, v2) == [[-3, -2], [-5]]
end

@testset "Is canonical" begin
    @test iscanonical(Int[]) == true
    @test iscanonical([-1]) == true
    @test iscanonical([-2, -2]) == true
    @test iscanonical([-3, -1, -1]) == true
    @test iscanonical([-1, -1, -2]) == false
    @test iscanonical([-1, -3]) == false
    @test iscanonical([-2, -2, -1, -1, -1]) == true
    @test iscanonical([-2, -2, 0, 0, 0]) == false
    @test iscanonical([1, -4, -4, -3]) == false
    @test iscanonical([-2, -1, -1]) == true
end

@testset "Canonical basis consistency" begin
    for n in 0:10
        @test iscanonical(-integer_partitions(n))
    end
end

@testset "Addition" begin
    v1 = Descendent([[-3, -2]], [1.0])
    v2 = Descendent([[-5], [-3, -2]], [1.0, 1.0])
    @test v1 + v2 == Descendent([[-3, -2], [-5]], [2.0, 1.0])
end

@testset "Once canonicalize" begin
    @test once_canonicalize([-1, -2]) == Descendent{Float64}([[-3], [-2, -1]], [1.0, 1.0])
    @test once_canonicalize([-1, -2, -1]) == Descendent{Float64}([[-3, -1], [-2, -1, -1]], [1.0, 1.0])
    @test once_canonicalize([-1, -2, -2]) == Descendent{Float64}([[-3, -2], [-5], [-2, -2, -1]], [2.0, 1.0, 1.0])
end

@testset "Descendent getindex" begin
    d = Descendent{Float64}([Int[], [1], [1, 1], [2], [3], [1, 2], [3, 2, 1]], [0.5, 1.5, 2.0, 2.1, 1.4, 1.3, 0.9])
    @test d[Int[]] == 0.5
    @test d[[1, 1]] == 2.0
    @test d[[2]] == 2.1
    @test d[[3, 2, 1]] == 0.9
end

reps_1 = vira_iter_solver(0.0, 1 / 2, 16)
reps_ψ = vira_iter_solver(1 / 2, 1 / 2, 16)
reps_σ = vira_iter_solver(1 / 16, 1 / 2, 16)


reps_ψ61 = Float64[7 1 0 0 4 0 0 6 0 0 0;
    0 6 3 4 0 0 0 0 0 0 0;
    0 0 5 0 8 1 0 9 0 0 0;
    0 0 0 5 0 3 9 0 3 0 0;
    0 0 0 0 0 4 0 9 4 0 0;
    0 0 0 0 0 0 4 0 6 16 0;
    0 0 0 0 0 0 0 0 0 3 36]


@testset "Level 2 rep" begin
    @test ConformalCalculator.vira_generator_level_col_solver(1 / 2, 1 / 2, [-2], 1, reps_ψ) == sparsevec([1], [3.0], 1)
    @test ConformalCalculator.vira_generator_level_col_solver(1 / 2, 1 / 2, [-1, -1], 1, reps_ψ) == sparsevec([1], [4.0], 1)
end

@testset "Virasoro iterate solver" begin
    @test reps_ψ[1] == [sparse([1.0;;])]
    @test reps_ψ[2] == [sparse([3.0 4.0]), sparse([2.25 3.0])]
    @test reps_ψ[3] == [sparse([4.0 1.0 0.0; 0.0 3.0 9.0]), sparse([5.0 6.25 15.0]), sparse([4.0 5.0 12.0])]
    @test reps_ψ[4] == [sparse([5.0 1.0 3.0 0.0 0.0; 0.0 4.0 6.0 4.0 0.0; 0.0 0.0 0.0 3.0 16.0]), sparse([6.0 0.0 12.5 3.0 0.0; 0.0 5.0 0.0 10.25 42.0]), sparse([7.0 10.0 15.0 20.0 72.0]), sparse([6.5 7.0 13.5 18.0 60.0])]
end
@testset "Level 6 rep" begin
    @test ConformalCalculator.vira_generator_level_col_solver(1 / 2, 1 / 2, [-2, -2, -2], 1, reps_ψ) == sparsevec([1, 3, 5], [6.0, 9.0, 9.0], 7)
    @test reps_ψ[6][1] == sparse(reps_ψ61)
end

@testset "Virasoro generator large matrix form" begin
    trunc_level = 4
    operator_level = 1
    M = ConformalCalculator._vira_block_embedding([reps_ψ[n][operator_level] for n in operator_level:trunc_level], operator_level, trunc_level)
    M_expect = spzeros(Float64, 12, 12)
    M_expect[1, 2] = 1.0
    M_expect[2, 3:4] = [3, 4]
    M_expect[3, 5:6] = [4, 1]
    M_expect[4, 6:7] = [3, 9]
    M_expect[5, 8:10] = [5, 1, 3]
    M_expect[6, 9:11] = [4, 6, 4]
    M_expect[7, 11:12] = [3, 16]
    @test M == M_expect
end

@testset "Virasoro representation property" begin
    L_rep_1 = vira_reps_total(reps_1)
    L_rep_ψ = vira_reps_total(reps_ψ)
    L_rep_σ = vira_reps_total(reps_σ)
    for L in [L_rep_1, L_rep_ψ, L_rep_ψ]
        for (m, Lm) in enumerate(L), (n, Ln) in enumerate(L)
            if m + n < length(L)
                @test Lm * Ln - Ln * Lm == (m - n) * L[m+n]
            end
        end
    end
end

expected_ranks_1 = [1, 0, 1, 1, 2, 2, 3, 3, 5, 5, 7, 8, 11, 12, 16, 18, 23]
expected_ranks_ψ = [1, 1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 9, 12, 14, 17, 20, 25]
expected_ranks_σ = [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32]
@testset "Gram matrix" begin
    gram_matrices_1 = gram_matrix(reps_1)
    gram_matrices_ψ = gram_matrix(reps_ψ)
    gram_matrices_σ = gram_matrix(reps_σ)
    @test map(rank, gram_matrices_1) == expected_ranks_1
    @test map(rank, gram_matrices_ψ)[1:13] == expected_ranks_ψ[1:13]
    @test map(rank, gram_matrices_σ)[1:13] == expected_ranks_σ[1:13]
end

@time ConformalCalculator.vira_generator_level_solver(0.0, 1 / 2, 11, 1, reps_1)