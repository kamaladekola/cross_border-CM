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
        zone_df[!, "Price"]  = results["λ"]["EOM"][z][end]

        for m in agents[:eom]
            if is_interconnector(m)
                z1, z2 = parse_agent_name(m) 
                if z == z1 || z == z2
                    zone_df[!, "G_$(m)"] = sign_for_zone(m, z) * results["g"][m][end]
                end
            else
                zone_m, _ = parse_agent_name(m)
                if zone_m == z
                    zone_df[!, "G_$(m)"] = results["g"][m][end]
                end
            end
        end
        CSV.write(joinpath(home_dir, "Results", "Scenario_$(scenario_overview_row["scen_number"])_EOM_Zone_$(z)_$(sens).csv"), zone_df; delim = ";")
    end
end



