# Save results
function save_results(mdict::Dict,EOM::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict,scenario_overview_row::DataFrameRow,sens) 
    # note that type of "sens" is not defined as a string stored in a dictionary is of type String31, whereas a "regular" string is of type String. Specifying one or the other may trow errors.
    vector_output = [scenario_overview_row["scen_number"]; sens; ADMM["n_iter"]; ADMM["walltime"];ADMM["Residuals"]["Primal"]["EOM"][end];ADMM["Residuals"]["Dual"]["EOM"][end]]
    CSV.write(joinpath(home_dir,string("overview_results.csv")), DataFrame(reshape(vector_output,1,:),:auto), delim=";",append=true);

    # EOM
    g_out = zeros(data["General"]["nTimesteps"],EOM["nAgents"])
    mm = 1
    for m in agents[:eom]
        g_out[:,mm] = results["g"][m][end]
        mm = mm+1
    end
    mat_output = [range(1,stop=data["General"]["nTimesteps"]) results[ "λ"]["EOM"][end] g_out -EOM["D"]]
    CSV.write(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_EOM_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Timestep";"Price";string.("G_",agents[:eom]);"Demand"]);
    
end
