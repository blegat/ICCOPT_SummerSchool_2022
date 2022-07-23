using JuMP
import SimpleConicADMM

function test_SOC()
    model = Model(SimpleConicADMM.Optimizer)
    set_optimizer_attributes(
        model,
        "ϵ_primal" => 1e-5,
        "ϵ_dual" => 1e-5,
        "ϵ_gap" => 1e-5,
    )
    @variable(model, x[1:3] in SecondOrderCone())
    @constraint(model, x[1] == 1)
    @objective(model, Max, x[2] + x[3])
    optimize!(model)
    display(solution_summary(model))
    return value.(x)
end

test_SOC()
