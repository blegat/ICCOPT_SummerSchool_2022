module TestSimpleConicADMM

using Test
using JuMP
import SimpleConicADMM

function test_runtests()
    optimizer = SimpleConicADMM.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iters"), 600)
    MOI.set(optimizer, MOI.Silent(), true) # comment this to enable output
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            optimizer,
        ),
        Float64,
    )
    config = MOI.Test.Config(
        atol = 1e-2,
        exclude = Any[
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.ConstraintName,
            MOI.VariableName,
            MOI.ObjectiveBound,
        ],
    )
    MOI.Test.runtests(
        model,
        config,
        exclude = String[
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
        ]
    )
    return
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

end  # module

TestSimpleConicADMM.runtests()
