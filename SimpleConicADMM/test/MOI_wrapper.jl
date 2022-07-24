module TestSimpleConicADMM

using Test
using JuMP
import SimpleConicADMM

function test_runtests()
    optimizer = SimpleConicADMM.Optimizer()
    #MOI.set(optimizer, MOI.Silent(), true) # uncomment this to suppress output
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
    MOI.Test.test_conic_SecondOrderCone_VectorAffineFunction(model, config)
    MOI.Test.test_conic_RotatedSecondOrderCone_VectorAffineFunction(model, config)
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
