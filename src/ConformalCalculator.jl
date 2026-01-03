module ConformalCalculator

using LinearAlgebra, SparseArrays
export O, HoloPoly, evaluation, derivative, Witt_action, Lie_action, deexponentialize, exp_Lie_action
export DescendentBasis, iscanonical, canonicalbasis, canonicalbasis_plain, once_canonicalize, Descendent, vira_iter_solver, gram_matrix
export CANONICALBASIS, CANONICALBASIS_PLAIN

include("./functionring.jl")
include("./deexponentialize.jl")
include("./conformalfamily.jl")
include("./canonicalbasis_buffer.jl")

end
