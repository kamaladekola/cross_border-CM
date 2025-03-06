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
        zone_df[!, "EOM_price"]  = results["λ"]["EOM"][z][end]
        zone_df[!, "CM_price"]  = fill(results["λ"]["CM"][z][end], nT)
        for m in agents[:eom]
            if m in agents[:IC]                                        # interconnector results
                zone_idx = findfirst(isequal(z), zones)
                zone_df[!, "$(m)"] = results["g"]["NetworkManager"][end][:, zone_idx]
            elseif m in agents[:Gen]                                                           # generator results
                zone_m, _ = parse_agent_name(m)
                if zone_m == z
                    zone_df[!, "$(m)"] = results["g"][m][end]
                    zone_df[!, "Capacity_$(m)"] = fill(results["y"][m][end], nT)
                    gen_name = get_agent_name(m)
                    zone_df[!, "new_capacity_$(m)"] = fill(results["y"][m][end] - data["Generators"][zone_m][gen_name]["C"], nT)
                end
            elseif m in agents[:Cons]
                cons_zone, _ = parse_agent_name(m)
                if cons_zone == z
                    zone_df[!, "$(m)"] = results["g"][m][end]
                    zone_df[!, "Inelastic_$(m)"] = results["Cons"]["inelastic_demand"][m][end] .* -1
                    zone_df[!, "Elastic_$(m)"] = results["Cons"]["elastic_demand"][m][end] .* -1
                    zone_df[!, "ENS_$(m)"] = results["Cons"]["ENS"][m][end]
                end
            end
        end
        CSV.write(joinpath(home_dir, "Results", "Scenario_$(scenario_overview_row["scen_number"])_EOM_Zone_$(z)_$(sens).csv"), zone_df; delim = ";")
    end
end



