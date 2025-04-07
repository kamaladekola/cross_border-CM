function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow, zones::Vector{String}, participation_matrix::DataFrame, derating_factor::DataFrame)

    # zone,_ = parse_agent_name(m)
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
    mod.ext[:parameters][:ρ_all] = data["ADMM"]["rho_EOM"] * ones(length(zones))

    # Parameters related to the crossborder electricity capacity markets
    mod.ext[:parameters][:λ_CM] = zeros(length(zones))
    mod.ext[:parameters][:cap_bar] = zeros(length(zones))
    mod.ext[:parameters][:ρ_CM] = data["ADMM"]["rho_CM"] * ones(length(zones))

    # ADMM parameters for interconnectors
    mod.ext[:parameters][:g_bar_all] = zeros(data["General"]["nTimesteps"], length(zones)) 
    mod.ext[:parameters][:λ_all] = zeros(data["General"]["nTimesteps"], length(zones))

    participation_switch = Dict{String, Dict{String, Float64}}()

    for m in eachrow(participation_matrix)
        agent_name = m[:agent]
        participation_switch[agent_name] = Dict(z => m[Symbol(z)] for z in zones)
    end

    data["General"]["participation_matrix"] = participation_switch
    mod.ext[:parameters][:participation_matrix] = data["General"]["participation_matrix"]

    # Derating factors
    derating_switch = Dict{String, Dict{String, Float64}}()

    for m in eachrow(derating_factor)
        agent_name = m[:agent]
        derating_switch[agent_name] = Dict(z => m[Symbol(z)] for z in zones)
    end

    data["General"]["derating_factor"] = derating_switch
    mod.ext[:parameters][:derating_factor] = data["General"]["derating_factor"]

    mod.ext[:parameters][:w] = ts[!, :weights][1:data["General"]["nTimesteps"]]

    return mod, agents
end