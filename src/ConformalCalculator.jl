module ConformalCalculator

using LinearAlgebra, SparseArrays, Combinatorics
export O, HoloPoly, evaluation, derivative, Witt_action, Lie_action, deexponentialize, exp_Lie_action
export iscanonical, once_canonicalize, Descendent, vira_iter_solver, vira_reps_total, gram_matrix

include("./functionring.jl")
include("./deexponentialize.jl")
include("./conformalfamily.jl")

end
