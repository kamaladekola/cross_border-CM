function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow)

    zone,_ = parse_agent_name(m)
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
    mod.ext[:parameters][:λ_EOM] = zeros(data["General"]["nTimesteps"])     # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["General"]["nTimesteps"])     # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["ADMM"]["rho_EOM"]                  # ADMM rho value

    # # Parameters related to the electricity CM
    mod.ext[:parameters][:λ_CM] = 0                                         # Price structure
    mod.ext[:parameters][:cap_bar] = 0                                      # ADMM penalty term
    mod.ext[:parameters][:ρ_CM] = data["ADMM"]["rho_CM"]                    # ADMM rho value

    mod.ext[:parameters][:CD] = data["CM"][zone]["capacity_target"]         # Get capacity target for zone Z
    mod.ext[:parameters][:WTP_CM] = data["CM"][zone]["price_target"]         # Willingness to pay for capacity in the CM
    mod.ext[:parameters][:CD_margin] = data["CM"][zone]["capacity_margin"] # Minimum willingness to pay for capacity in the CM

    

    return mod, agents
end