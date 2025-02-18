function ADMM_subroutine!(m::String,results::Dict,ADMM::Dict,EOM::Dict, CM:: Dict, mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
zone, _ = parse_agent_name(m)

# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents_z"]+1)*last(ADMM["Imbalances"]["EOM"][zone])
    mod.ext[:parameters][:λ_EOM] = last(results["λ"]["EOM"][zone])
    mod.ext[:parameters][:ρ_EOM] = last(ADMM["ρ"]["EOM"][zone])

    if m in agents[:cm]
        mod.ext[:parameters][:cap_bar] = results["cap_cm"][m][end] - 1/(CM["nAgents_z"]+1)*last(ADMM["Imbalances"]["CM"][zone])
        mod.ext[:parameters][:λ_CM] = last(results["λ"]["CM"][zone])
        mod.ext[:parameters][:ρ_CM] = last(ADMM["ρ"]["CM"][zone])
    end 
end

# Solve agents decision problems:
if m in agents[:Gen]
    @timeit TO_local "Solve generator problems" begin
        solve_generator_agent!(mod)  
    end
elseif m in agents[:Cons]
    @timeit TO_local "Solve consumer problems" begin
        solve_consumer_agent!(mod)  
    end

elseif m in agents[:IC]
    z1, z2 = parse_agent_name(m)
    mod.ext[:parameters][:λ1] = last(results["λ"]["EOM"][z1])
    mod.ext[:parameters][:λ2] = last(results["λ"]["EOM"][z2])
    @timeit TO_local "Solve interconnector problems" begin
        solve_interconnector_agent!(mod)
    end
end
# Query results
@timeit TO_local "Query results" begin
    push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))                                                # generation
                                                  
    if m in agents[:Gen]
        push!(results["y"][m], value(mod.ext[:variables][:y]))                                                      # installed capacity
    end

    if m in agents[:cm]
        push!(results["cap_cm"][m], value.(mod.ext[:variables][:cap_cm]))                                           # capacity offered in capacity markets
    end

    if m in agents[:Cons]
        push!(results["Cons"]["inelastic_demand"][m], collect(value.(mod.ext[:variables][:g_VOLL])))                        # inelastic demand
        push!(results["Cons"]["elastic_demand"][m], collect(value.(mod.ext[:variables][:g_ela])))                           # elastic demand
        push!(results["Cons"]["ENS"][m], collect(value.(mod.ext[:variables][:ens])))                                        # unserved energy
    end
    
end

# Merge local TO with TO:
merge!(TO,TO_local)
end