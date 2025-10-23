function save_results(mdict::Dict, EOM::Dict, ADMM::Dict, results::Dict, data::Dict, agents::Dict, scenario_overview_row::DataFrameRow,sens, zones::Vector{String})
    # note that type of "sens" is not defined as a string stored in a dictionary is of type String31, whereas a "regular" string is of type String. Specifying one or the other may throw errors.
    primal_residuals = [ADMM["Residuals"]["Primal"]["EOM"][z][end] for z in zones]
    dual_residuals   = [ADMM["Residuals"]["Dual"]["EOM"][z][end]   for z in zones]

    vector_output = vcat(scenario_overview_row["scen_number"], sens, ADMM["n_iter"], ADMM["walltime"], primal_residuals..., dual_residuals...)

    overview_header = vcat(["scen_number", "sensitivity", "n_iter", "walltime"], ["PrimalResidual_EOM_$z" for z in zones], ["DualResidual_EOM_$z"   for z in zones])

    overview_path = joinpath(home_dir, "overview_results.csv")
    df_overview   = DataFrame(reshape(vector_output, 1, :), :auto)
    rename!(df_overview, overview_header)

    if isfile(overview_path)
        CSV.write(overview_path, df_overview; delim = ";", append = true)
    else
        CSV.write(overview_path, df_overview; delim = ";")
    end

    # Results for each zone
    nT = data["General"]["nTimesteps"]
    timesteps = 1:nT

    for z in zones
        zone_df = DataFrame(Timestep = timesteps)
        zone_df[!, "EOM_price"] = results["λ"]["EOM"][z][end]
        zone_df[!, "CM_price"] = fill(results["λ"]["CM"][z][end], nT)
        
        zone_idx = findfirst(isequal(z), zones)
        
        # Add total capacity offered/demanded in zone
        total_cap_offered = sum(results["cap_cm"][m][end][zone_idx] for m in agents[:cm] if m in agents[:Gen])
        total_local_cap_offered = sum(sum(results["cap_cm"][m][end]) for m in agents[:cm_Z][z] if m in agents[:Gen])
        total_cap_demanded = sum(results["cap_cm"][m][end][zone_idx] for m in agents[:cm] if m in agents[:Cons])
        zone_df[!, "TotalCapOffer"] = fill(total_cap_offered, nT)
        zone_df[!, "TotalCapDemanded"] = fill(total_cap_demanded, nT)
        zone_df[!, "TotalLocalCapSupply"] = fill(total_local_cap_offered, nT)
        
        for m in agents[:all]
            if m in agents[:IC]
                zone_idx = findfirst(isequal(z), zones)
                zone_df[!, "$(m)"] = results["g"]["NetworkManager"][end][:, zone_idx]
            elseif m in agents[:CIC]
                zone_idx = findfirst(isequal(z), zones)
                zone_df[!, "$(m)"] = fill(results["cap_cm"][m][end][zone_idx], nT)
            elseif m in agents[:Gen]
                zone_m, _ = parse_agent_name(m)
                if zone_m == z
                    zone_df[!, "$(m)"] = results["g"][m][end]
                    zone_df[!, "new_capacity_$(m)"] = fill(results["y"][m][end], nT)
                    gen_name = get_agent_name(m)
                    zone_df[!, "Capacity_$(m)"] = fill(results["y"][m][end] + data["Generators"][zone_m][gen_name]["C"], nT)
                    
                    if m in agents[:cm]
                        for (tgt_idx, target_z) in enumerate(zones)
                            # zone_df[!, "Cap_to_$(target_z)_$(m)"] = fill(results["cap_cm"][m][end][tgt_idx], nT)
                            # Sum capacity offered FROM current zone z TO target zone
                            local_cap_to_target = sum(results["cap_cm"][m][end][tgt_idx] for m in agents[:cm_Z][z] if m in agents[:Gen])
                            # zone_df[!, "Offered_Cap_from_$(z)_to_$(target_z)"] = fill(local_cap_to_target, nT)
                        end
                    end
                end
            elseif m in agents[:Cons]
                cons_zone, _ = parse_agent_name(m)
                if cons_zone == z
                    # Existing consumer results
                    zone_df[!, "$(m)"] = results["g"][m][end]
                    zone_df[!, "Inelastic_$(m)"] = results["Cons"]["inelastic_demand"][m][end] .* -1
                    zone_df[!, "Elastic_$(m)"] = results["Cons"]["elastic_demand"][m][end] .* -1
                    zone_df[!, "ENS_$(m)"] = results["Cons"]["ENS"][m][end]
                    
                    if m in agents[:cm]
                        zone_df[!, "CapDemand_$(m)"] = fill(results["cap_cm"][m][end][zone_idx], nT)
                    end
                end
            end
        end
        
        CSV.write(joinpath(home_dir, "Results", "Scenario_$(scenario_overview_row["scen_number"])_EOM_Zone_$(z)_$(sens).csv"), zone_df; delim = ";")
    end
end



