module SimpleConicADMM

using LinearAlgebra, SparseArrays

import MathOptInterface
const MOI = MathOptInterface

include("solver.jl")
include("MOI_wrapper.jl")

end # module
