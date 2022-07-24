const SUPPORTED_CONES = Union{
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
}

MOI.Utilities.@product_of_sets(
    Cones,
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
)

const OptimizerCache{T} = MOI.Utilities.GenericModel{
    T,
    MOI.Utilities.ObjectiveContainer{T},
    MOI.Utilities.VariablesContainer{T},
    MOI.Utilities.MatrixOfConstraints{
        T,
        MOI.Utilities.MutableSparseMatrixCSC{
            T,
            Int,
            MOI.Utilities.OneBasedIndexing,
        },
        Vector{T},
        Cones{T},
    },
}

const DEFAULT_OPTIONS = Dict{String,Any}(
    "max_iters" => 100,
    "ϵ_primal" => 1e-4,
    "ϵ_dual" => 1e-4,
    "ϵ_gap" => 1e-4,
    "ϵ_unbounded" => 1e-7,
    "ϵ_infeasible" => 1e-7,
)

"""
    Optimizer()

Create a new ECOS optimizer.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    cones::Union{Nothing,Cones{T}}
    data::Union{Nothing,Data{T}}
    sol::Union{Nothing,Solution{T}}
    cache::Union{Nothing,Cache{T}}
    max_sense::Bool
    objective_constant::T
    silent::Bool
    options::Dict{String,Any}

    function Optimizer{T}() where T
        return new{T}(nothing, nothing, nothing, nothing, false, zero(T), false, copy(DEFAULT_OPTIONS))
    end
end
Optimizer() = Optimizer{Float64}()

function MOI.default_cache(::Optimizer{T}, ::Type{T}) where {T}
    return MOI.Utilities.UniversalFallback(OptimizerCache{T}())
end

function MOI.is_empty(optimizer::Optimizer)
    return optimizer.cones === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.cones = nothing
    optimizer.data = nothing
    optimizer.sol = nothing
    optimizer.cache = nothing
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "SimpleConicADMM"

MOI.get(::Optimizer, ::MOI.SolverVersion) = "v0.0.1"

# MOI.RawOptimizerAttribute

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return hasfield(settings, Symbol(param.name))
end

function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
    optimizer.options[param.name] = value
    return
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    if !haskey(optimizer.options, param.name)
        throw(MOI.UnsupportedAttribute(param))
    end
    return optimizer.options[param.name]
end

# MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

# MOI.supports

function MOI.supports(
    ::Optimizer{T},
    ::Union{
        MOI.ObjectiveSense,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    },
) where T
    return true
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.VectorAffineFunction{T}},
    ::Type{<:SUPPORTED_CONES},
) where T
    return true
end

function MOI.optimize!(optimizer::Optimizer{T}) where T
    data = optimizer.data
    sol = optimizer.sol
    options = optimizer.options
    max_iters = options["max_iters"]
    primal_feasible = false
    dual_feasible = false
    norm_b = norm(data.b)
    norm_c = norm(data.c)
    n = length(data.c)
    start = time()
    sol.solve_time = @elapsed for i in 0:max_iters
        x = primal(sol, n)
        y = dual(sol, n)
        p = (data.A * x + slack(sol, n)) - data.b
        d = data.A' * y + data.c
        cx = dot(data.c, x)
        by = dot(data.b, y)
        primal_feasibility = norm(p) / (one(T) + norm_b)
        primal_feasible = primal_feasibility ≤ options["ϵ_primal"]
        dual_feasibility = norm(d) / (one(T) + norm_c)
        dual_feasible = dual_feasibility ≤ options["ϵ_dual"]
        g = cx + by
        rel_gap = abs(g) / (1 + abs(cx) + abs(by))
        if !optimizer.silent
            print_info(i, primal_feasibility, cx, dual_feasibility, by, rel_gap, start)
        end
        if primal_feasible && dual_feasible && rel_gap ≤ options["ϵ_gap"]
            sol.raw_status = "Solved to optimality"
            sol.termination_status = MOI.OPTIMAL
            sol.primal_status = MOI.FEASIBLE_POINT
            sol.dual_status = MOI.FEASIBLE_POINT
            return
        end
        ux = unscaled_primal(sol, n)
        vs = unscaled_slack(sol, n)
        cux = -dot(data.c, ux)
        if norm(data.A * ux + vs) < (cux / norm_c) * options["ϵ_unbounded"]
            sol.raw_status = "Detected to be dual infeasible"
            sol.termination_status = MOI.DUAL_INFEASIBLE
            sol.primal_status = MOI.INFEASIBILITY_CERTIFICATE
            sol.dual_status = MOI.INFEASIBLE_POINT
            sol.uxy /= cux
            sol.vrs /= cux
            return
        end
        uy = unscaled_dual(sol, n)
        buy = -dot(data.b, uy)
        if norm(data.A' * uy) < (buy / norm_b) * options["ϵ_infeasible"]
            sol.raw_status = "Detected to be primal infeasible"
            sol.termination_status = MOI.INFEASIBLE
            sol.primal_status = MOI.INFEASIBLE_POINT
            sol.dual_status = MOI.INFEASIBILITY_CERTIFICATE
            sol.uxy /= buy
            sol.vrs /= buy
            return
        end
        sol.iter = i
        if i >= max_iters
            break
        end
        iterate!(sol, optimizer.data, optimizer.cache, optimizer.cones, size(data.A, 2))
    end
    sol.raw_status = "Maximum number of iterations ($max_iters) reached"
    sol.termination_status = MOI.ITERATION_LIMIT
    sol.primal_status = primal_feasible ? MOI.FEASIBLE_POINT : MOI.INFEASIBLE_POINT
    sol.dual_status = dual_feasible ? MOI.FEASIBLE_POINT : MOI.INFEASIBLE_POINT
    return
end

function MOI.copy_to(dest::Optimizer{T}, src::OptimizerCache{T}) where {T}
    MOI.empty!(dest)
    Ab = src.constraints
    A = -convert(SparseMatrixCSC{T,Int}, Ab.coefficients)
    dest.max_sense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    obj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
    dest.objective_constant = MOI.constant(obj)
    c0 = zeros(A.n)
    for term in obj.terms
        c0[term.variable.value] += term.coefficient
    end
    dest.cones = deepcopy(src.constraints.sets)
    dest.data = Data(
        A,
        Ab.constants,
        dest.max_sense ? -c0 : c0,
    )
    dest.sol, dest.cache = setup(dest.data)
    return MOI.Utilities.identity_index_map(src)
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache},
)
    MOI.Utilities.throw_unsupported(src)
    return MOI.copy_to(dest, src.model)
end

function MOI.copy_to(dest::Optimizer{T}, src::MOI.ModelLike) where {T}
    cache = OptimizerCache{T}()
    index_map = MOI.copy_to(cache, src)
    MOI.copy_to(dest, cache)
    return index_map
end

MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec) = optimizer.sol.solve_time

function MOI.get(optimizer::Optimizer, ::MOI.BarrierIterations)
    return optimizer.sol.iter
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    if optimizer.sol === nothing
        return RAW_OPTIMIZE_NOT_CALLED
    else
        return optimizer.sol.raw_status
    end
end

# Implements getter for result value and statuses
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    return isnothing(optimizer.sol) ? MOI.OPTIMIZE_NOT_CALLED : optimizer.sol.termination_status
end

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    c = optimizer.data.c
    n = length(c)
    value = dot(c, unscaled_primal(optimizer.sol, n)) / optimizer.sol.uτ
    if optimizer.max_sense
        value = -value
    end
    if !MOI.Utilities.is_ray(MOI.get(optimizer, MOI.PrimalStatus()))
        value += optimizer.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    b = optimizer.data.b
    n = length(optimizer.data.c)
    value = dot(b, unscaled_dual(optimizer.sol, n)) / optimizer.sol.uτ
    if !optimizer.max_sense
        value = -value
    end
    if !MOI.Utilities.is_ray(MOI.get(optimizer, MOI.DualStatus()))
        value += optimizer.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return optimizer.sol.primal_status
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return primal(optimizer.sol, length(optimizer.data.c))[vi.value]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return slack(optimizer.sol, length(optimizer.data.c))[MOI.Utilities.rows(optimizer.cones, ci)]
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return optimizer.sol.dual_status
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return dual(optimizer.sol, length(optimizer.data.c))[MOI.Utilities.rows(optimizer.cones, ci)]
end

function MOI.get(optimizer::Optimizer, ::MOI.ResultCount)
    if isnothing(optimizer.sol)
        return 0
    else
        return 1
    end
end
