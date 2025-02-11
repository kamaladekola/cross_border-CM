function ADMM_subroutine!(m::String,results::Dict,ADMM::Dict,EOM::Dict,mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
zone, _ = parse_agent_name(m)

# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents_z"]+1)*last(ADMM["Imbalances"]["EOM"][zone])
    mod.ext[:parameters][:λ_EOM] = last(results["λ"]["EOM"][zone])
    mod.ext[:parameters][:ρ_EOM] = last(ADMM["ρ"]["EOM"][zone])
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
    push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
end

# Merge local TO with TO:
merge!(TO,TO_local)
end