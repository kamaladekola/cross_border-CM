function ADMM_subroutine!(m::String,results::Dict,ADMM::Dict,EOM::Dict,mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][end]
    mod.ext[:parameters][:λ_EOM] = results["λ"]["EOM"][end] 
    mod.ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
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
end

# Query results
@timeit TO_local "Query results" begin
    push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
end

# Merge local TO with TO:
merge!(TO,TO_local)
end