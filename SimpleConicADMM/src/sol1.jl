struct Data{T}
    A::SparseMatrixCSC{T,Int}
    b::Vector{T}
    c::Vector{T}
end


const RAW_OPTIMIZE_NOT_CALLED = "Optimize not called"
mutable struct Solution{T}
    uxy::Vector{T}
    uτ::T
    vrs::Vector{T}
    vκ::T
    raw_status::String
    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    solve_time::Float64
    iter::Int
    function Solution{T}(k) where T
        return new{T}(
            zeros(T, k),
            one(T),
            zeros(T, k),
            one(T),
            RAW_OPTIMIZE_NOT_CALLED,
            MOI.OPTIMIZE_NOT_CALLED,
            MOI.UNKNOWN_RESULT_STATUS,
            MOI.UNKNOWN_RESULT_STATUS,
            0.0,
            0,
        )
    end
end

function unscaled_primal(sol::Solution, n)
    return view(sol.uxy, 1:n)
end
function primal(sol::Solution, n)
    return unscaled_primal(sol, n) / sol.uτ
end

function unscaled_dual(sol::Solution, n)
    return view(sol.uxy, (n + 1):length(sol.uxy))
end
function dual(sol::Solution, n)
    return unscaled_dual(sol, n) / sol.uτ
end

function unscaled_slack(sol::Solution, n)
    return view(sol.vrs, (n + 1):length(sol.vrs))
end
function slack(sol::Solution, n)
    return unscaled_slack(sol, n) / sol.uτ
end

import SuiteSparse
struct Cache{T}
    IQ::SparseMatrixCSC{T,Int}
end

function setup(data::Data{T}) where {T}
    m, n = size(data.A)
    A = data.A
    b = data.b
    c = data.c
    IM = [I A'; -A I]
    h = [c; b]
    IQ = [
        IM h
        -h' 1
    ]
    # TODO Precomputes quantitiees here
    solution = Solution{T}(n + m)
    cache = Cache{T}(
        IQ
    )
    return solution, cache
end

function project_affine(data, cache, wxy, wτ, n)
    w = [wxy; wτ]
    u = cache.IQ \ w
    uxy = u[1:end-1]
    uτ = u[end]
    return uxy, uτ
end

function _set(cones, ci::MOI.ConstraintIndex{F,S}) where {F,S}
    d = length(MOI.Utilities.rows(cones, ci))
    return MOI.Utilities.set_with_dimension(S, d)
end

import MathOptSetDistances
function _project_cones(cones, u, w, cis, n)
    for ci in cis
        cone = _set(cones, ci)
        dual = MOI.dual_set(cone)
        rows = MOI.Utilities.rows(cones, ci)
        u[n .+ rows] = MathOptSetDistances.projection_on_set(
            MathOptSetDistances.DefaultDistance(),
            w[n .+ rows],
            dual,
        )
    end
end

function project_cones(cones, w, n)
    u = copy(w)
    for (F, S) in MOI.get(cones, MOI.ListOfConstraintTypesPresent())
        cis = MOI.get(cones, MOI.ListOfConstraintIndices{F,S}())
        _project_cones(cones, u, w, cis, n)
    end
    return u
end

function iterate!(sol::Solution{T}, data, cache, cones, n) where {T}
    ũxy, ũτ = project_affine(data, cache, sol.uxy + sol.vrs, sol.uτ + sol.vκ, n)
    sol.uxy = project_cones(cones, ũxy - sol.vrs, n)
    sol.uτ = max(zero(T), ũτ - sol.vκ)
    sol.vrs = sol.vrs - ũxy + sol.uxy
    sol.vκ = sol.vκ - ũτ + sol.uτ
    return
end

using Printf

# print objective gap information for iterative
function print_info(i, primal_feasibility, cx, dual_feasibility, by, rel_gap, start)
    if iszero(i)
        @printf "\n%-5s | %-14s | %-14s | %-14s | %-14s | %-11s | %-11s\n" "Iter." "Primal Feas." "Primal Obj." "Dual Feas." "Dual Obj." "Rel. gap" "Time (s)"
    end
    @printf "%5d | %+14.6e | %+14.6e | %+14.6e | %+14.6e | %11.3e | %11.3e\n" i primal_feasibility cx dual_feasibility -by rel_gap (time() - start)
    flush(stdout)
    flush(stderr)
    return
end
