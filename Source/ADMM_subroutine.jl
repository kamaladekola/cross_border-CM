function ADMM_subroutine!(m::String, results::Dict, ADMM::Dict, EOM::Dict, CM::Dict, mod::Model, agents::Dict, TO::TimerOutput, zones::Vector{String})
    TO_local = TimerOutput()

    if m == "NetworkManager"
        # Update network manager-specific parameters
        @timeit TO_local "Compute NetworkManager penalty terms" begin
            n_Rows = size(results["g"][m][end], 1)
            n_Zones = length(zones)
            mod.ext[:parameters][:g_bar_all] = Matrix{Float64}(undef, n_Rows, n_Zones)
            for (zone_idx, z) in enumerate(zones)
                mod.ext[:parameters][:g_bar_all][:, zone_idx] = results["g"][m][end][:, zone_idx] .- (1/(EOM["nAgents_z"][z]+1)) * last(ADMM["Imbalances"]["EOM"][z])
            end
            mod.ext[:parameters][:λ_all] = hcat([last(results["λ"]["EOM"][z]) for z in zones]...)
            mod.ext[:parameters][:ρ_all] = [last(ADMM["ρ"]["EOM"][z]) for z in zones]
        end

        @timeit TO_local "Solve network manager problems" begin
            solve_interconnector_agent!(mod)
        end
    else
        # i.e Gen and Cons agents
        zone, _ = parse_agent_name(m)
        @timeit TO_local "Compute ADMM penalty terms" begin
            mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents_z"][zone]+1) * last(ADMM["Imbalances"]["EOM"][zone])
            mod.ext[:parameters][:λ_EOM] = last(results["λ"]["EOM"][zone])
            mod.ext[:parameters][:ρ_EOM] = last(ADMM["ρ"]["EOM"][zone])
    
            # Update CM penalty terms
            if m in agents[:cm]
                for (zone_idx, z) in enumerate(zones)
                    mod.ext[:parameters][:cap_bar][zone_idx] = results["cap_cm"][m][end][zone_idx] - (1/(CM["nAgents_z"][z]+1)) * last(ADMM["Imbalances"]["CM"][z])
                end
                mod.ext[:parameters][:λ_CM] = hcat([last(results["λ"]["CM"][zone]) for zone in zones]...)
                mod.ext[:parameters][:ρ_CM] = [last(ADMM["ρ"]["CM"][zone]) for zone in zones]
            end 
        end
    
        if m in agents[:Gen]
            @timeit TO_local "Solve generator problems" begin
                solve_generator_agent!(mod, m, zones)
            end
        elseif m in agents[:Cons]
            @timeit TO_local "Solve consumer problems" begin
                solve_consumer_agent!(mod, m, zones)
            end
        end
    end

    # Query results block (update accordingly per agent type)
    @timeit TO_local "Query results" begin
        if m in agents[:Gen]
            push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
            push!(results["y"][m], value(mod.ext[:variables][:y]))
            # Update capacity market results if this agent participates in CM
            if m in agents[:cm]
                push!(results["cap_cm"][m], collect(value.(mod.ext[:variables][:cap_cm])))
            end
        elseif m in agents[:Cons]
            push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
            push!(results["Cons"]["inelastic_demand"][m], collect(value.(mod.ext[:variables][:g_VOLL])))
            push!(results["Cons"]["elastic_demand"][m], collect(value.(mod.ext[:variables][:g_ela])))
            push!(results["Cons"]["ENS"][m], collect(value.(mod.ext[:variables][:ens])))
            # Update capacity market results if this agent participates in CM
            if m in agents[:cm]
                push!(results["cap_cm"][m], collect(value.(mod.ext[:variables][:cap_cm])))
            end
        elseif m == "NetworkManager"
            push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
        end
    end

    merge!(TO, TO_local)
end
