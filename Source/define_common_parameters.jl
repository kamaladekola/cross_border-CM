function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow, zones::Vector{String})

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
    mod.ext[:sets][:JZ] = 1:length(zones)
    mod.ext[:sets][:JL] = 1:data["Network"]["nLines"]
  
    # Parameters related to the EOM
    mod.ext[:parameters][:λ_EOM] = zeros(data["General"]["nTimesteps"])     # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["General"]["nTimesteps"])     # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["ADMM"]["rho_EOM"]                  # ADMM rho value
    mod.ext[:parameters][:ρ_all] = ones(data["General"]["nTimesteps"])

    # # Parameters related to the electricity CM
    mod.ext[:parameters][:λ_CM] = 0                                         # Price structure
    mod.ext[:parameters][:cap_bar] = 0                                      # ADMM penalty term
    mod.ext[:parameters][:ρ_CM] = data["ADMM"]["rho_CM"]                    # ADMM rho value

    # ADMM parameters for interconnectors
    mod.ext[:parameters][:g_bar_all] = zeros(data["General"]["nTimesteps"], length(zones)) 
    mod.ext[:parameters][:λ_all] = zeros(data["General"]["nTimesteps"], length(zones))
    return mod, agents
end