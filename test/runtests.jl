using Test
using FunctionRing

println("---------")
println("O")
println("---------")

setup = [
    O{Int}(2),
    O{Int}(3),
    O{Int}(Inf),
    O{Rational}(2 // 3),
    O{Rational}(Inf),
    O{Rational}(3.5)
]

relations_true = [
    2 < O{Int}(3),
    4 > O{Int}(3),
    2 < O{Int}(3.5),
    2 < O{Int}(Inf),
    3 < O{Rational}(3.5),
    4 > O{Rational}(3.4),
    3 < O{Rational}(Inf),
    -2 < O{Int}(3),
]

relations_false = [
    2 < O{Int}(1),
    2 > O{Int}(3.5),
    2 > O{Int}(Inf),
    3 > O{Rational}(3.5),
]

equations = [
    O{Int}(2) + O{Int}(3),
    O{Rational}(2) + O{Rational}(2 // 3),
    3 * O{Rational}(4),
    3 // 2 * O{Rational}(4),
    4.5 * O{Rational}(4),
    [2, 3, 4] * O{Int}(2),
    [2, 3, 4] * O{Rational}(2.5),
    [3, 4, 2] * O{Rational}(2.5)
]

answers = [
    O{Int}(2),
    O{Rational}(2 // 3),
    O{Rational}(7),
    O{Rational}(4 + 3 // 2),
    O{Rational}(4 + 4.5),
    O{Int}(4),
    O{Rational}(2.5 + 2),
    O{Rational}(2.5 + 2)
]

@testset "Setup" begin
    for t in setup
        t
    end

    @test true
end

@testset "True relations" begin
    for rel in relations_true
        rel == true
    end
end

@testset "False relations" begin
    for rel in relations_false
        rel == false
    end
end

@testset "Identities" begin
    for (eq, ans) in zip(equations, answers)
        @test eq == ans
    end
end

println("-------------------------------------------------------")
println("Function construction: Setup, reordering and truncation")
println("-------------------------------------------------------")

setup_fun = [
    HoloPoly{Int,Float64}([0, 1, 3, 4], [2.4, 3.3, 4.0, 5.6]),
    HoloPoly{Int,Float64}([1, 3, 4], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0]),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(Inf)),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3)),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2)),
    HoloPoly{Int,ComplexF64}([0, 1, 3, 4], [2.4 + im, 3.3, 4.0, 5.6]),
    HoloPoly{Int,ComplexF64}([1, 3, 4], [2.4, 3.3, -4.0 + 0.3 * im]),
    HoloPoly{Int,ComplexF64}([3, 4, 1], [2.4, 3.3 + 0.2 * im, -4.0]),
    HoloPoly{Int,ComplexF64}([-2, 4, 1], [2.4, 3.3 + 0.2 * im, -4.0]),
    HoloPoly{Rational{Int64},ComplexF64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3 + 5 * im, -4.0]),
]

power_reordering = [
    HoloPoly{Int,Float64}(Int[1], Float64[1e-15]),
    HoloPoly{Int,Float64}(Int[1, 3, 5, 2], Float64[1e-15, 2*1e-15, 0.4*1e-16, -0.3*1e-20]),
    HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([-3, 4, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([-3, 4, 1], [2.4, 3.3, -4.0], O{Int}(3)),
    HoloPoly{Int,Float64}([-3, 4, 1], [2.4, 3.3, -4.0], O{Int}(-4)),
    HoloPoly{Int,ComplexF64}([-3, 4, 1], [2.4, 3.3, -4.0 + 2im]),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2)),
]

power_reordering_answer = [
    HoloPoly{Int,Float64}(Int[], Float64[]),
    HoloPoly{Int,Float64}(Int[], Float64[]),
    HoloPoly{Int,Float64}([1, 3, 4], [-4.0, 2.4, 3.3]),
    HoloPoly{Int,Float64}([-3, 1, 4], [2.4, -4.0, 3.3]),
    HoloPoly{Int,Float64}([-3, 1], [2.4, -4.0], O{Int}(3)),
    HoloPoly{Int,Float64}(Int[], Float64[], O{Int}(-4)),
    HoloPoly{Int,ComplexF64}([-3, 1, 4], [2.4, -4.0 + 2im, 3.3]),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4], [2.4, 3.3], O{Rational{Int64}}(3 // 2)),
]

@testset "Setup" begin
    for t in setup_fun
        t
    end

    @test true
end

@testset "Reordering and truncation" begin
    for (t, ans) in zip(power_reordering, power_reordering_answer)
        @test t == ans
    end
end

println("------------------------------")
println("Function calculation: addition")
println("------------------------------")

addition = [
    HoloPoly{Int,Float64}(Int[1], Float64[1e-15]) + HoloPoly{Int,Float64}(Int[1, 3, 5, 2], Float64[1e-15, 2*1e-15, 0.4*1e-16, -0.3*1e-20]),
    HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]) + HoloPoly{Int,Float64}([-3, -2, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([1, 2, 3], [2.4, 3.3, -4.0]) + HoloPoly{Int,Float64}([-3, -2, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([1, 2, 3], [2.4, 3.3, -4.0]) - HoloPoly{Int,Float64}([-3, -2, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,ComplexF64}([-3, 4, 1], [2.4, 3.3, -4.0 + 2im]) + HoloPoly{Int,ComplexF64}([-1, 2, 1], [2.4, 3.3, -2.0 + 2im]),
    HoloPoly{Int,ComplexF64}([-1, 0, 1], [2.4, 3.3, -2.0 + 2im]) + HoloPoly{Int,ComplexF64}([0, 2, 3], [2.4, 3.3, -2.0 + 2im], O{Int}(2)),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2)) + HoloPoly{Rational{Int64},Float64}([1 // 16, 3 // 4, 3 // 2], [2.4, 3.3, -4.0], O{Rational{Int64}}(4 // 2)),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2)) + HoloPoly{Rational{Int64},Float64}([1 // 1], [1.0], O{Rational{Int64}}(1 // 3)),
    HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2)) - HoloPoly{Rational{Int64},Float64}([1 // 2, 3 // 4, 4 // 1], [2.4, 3.3, -4.0], O{Rational{Int64}}(3 // 2))
]

sum_result = [
    HoloPoly{Int,Float64}(Int[], Float64[]),
    HoloPoly{Int,Float64}([-3, -2, 1, 3, 4], [2.4, 3.3, -8.0, 2.4, 3.3]),
    HoloPoly{Int,Float64}([-3, -2, 1, 2, 3], [2.4, 3.3, 2.4 - 4.0, 3.3, -4.0]),
    HoloPoly{Int,Float64}([-3, -2, 1, 2, 3], [-2.4, -3.3, 2.4 + 4.0, 3.3, -4.0]),
    HoloPoly{Int,ComplexF64}([-3, -1, 1, 2, 4], [2.4, 2.4, -6.0 + 4im, 3.3, 3.3]),
    HoloPoly{Int,ComplexF64}([-1, 0, 1], [2.4, 5.7, -2.0 + 2im], O{Int}(2)),
    HoloPoly{Rational{Int64},Float64}([1 // 16, 1 // 2, 3 // 4], [2.4, 2.4, 6.6], O{Rational{Int64}}(3 // 2)),
    HoloPoly{Rational{Int64},Float64}(Rational{Int64}[], Float64[], O{Rational{Int64}}(1 // 3)),
    HoloPoly{Rational{Int64},Float64}(Rational{Int64}[], Float64[], O{Rational{Int64}}(3 // 2))
]

@testset "Addition" begin
    for (add, result) in zip(addition, sum_result)
        @test add == result
    end
end

println("-------------------------------------------")
println("Function calculation: number multiplication")
println("-------------------------------------------")

number_multiplication = [
    2.5 * HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    1.0 * HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    -4.0 * HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    1.5 * HoloPoly{Int,Float64}([0], [0.0], O{Int}(2)),
    (1.5 + im) * HoloPoly{Rational{Int64},ComplexF64}([1 // 2, 1 // 16, 3 // 2], [0.5 + 0.5im, 1.5im, 2.0], O{Rational{Int64}}(5 // 4))
]

number_multiplication_ans = [
    HoloPoly{Int,Float64}([3, 4, 1], [2.5 * 2.4, 2.5 * 3.3, -2.5 * 4.0]),
    HoloPoly{Int,Float64}([3, 4, 1], [2.4, 3.3, -4.0]),
    HoloPoly{Int,Float64}([3, 4, 1], [-4 * 2.4, -4 * 3.3, 4 * 4.0]),
    HoloPoly{Int,Float64}(Int[], Float64[], O{Int}(2)),
    HoloPoly{Rational{Int64},ComplexF64}([1 // 2, 1 // 16], [(0.5 + 0.5im) * (1.5 + im), 1.5im * (1.5 + im)], O{Rational{Int64}}(5 // 4))
]

@testset "Number multiplication" begin
    for (number_mult, result) in zip(number_multiplication, number_multiplication_ans)
        @test number_mult == result
    end
end

println("---------------------------------------------")
println("Function calculation: function multiplication")
println("---------------------------------------------")

