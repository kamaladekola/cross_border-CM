function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow)
    # Solver settings
    # Define dictonaries for sets, parameters, timeseries, variables, constraints & expressions
    mod.ext[:sets] = Dict()
    mod.ext[:parameters] = Dict()
    mod.ext[:timeseries] = Dict()
    mod.ext[:variables] = Dict()
    mod.ext[:constraints] = Dict()
    mod.ext[:expressions] = Dict()

    # Sets
    mod.ext[:sets][:JH] = 1:data["General"]["nTimesteps"]
  
    # Parameters related to the EOM
    mod.ext[:parameters][:λ_EOM] = zeros(data["General"]["nTimesteps"])   # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["General"]["nTimesteps"])   # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["ADMM"]["rho_EOM"]                # ADMM rho value 

    return mod, agents
end