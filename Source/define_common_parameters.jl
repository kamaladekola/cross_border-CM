function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow, zones::Vector{String}, ptdf::DataFrame, nodal_ptdf::DataFrame, lines::DataFrame, participation_matrix::DataFrame, derating_factor::DataFrame)

    # zone, _ = parse_agent_name(m)


    mod.ext[:sets] = Dict()
    mod.ext[:parameters] = Dict()
    mod.ext[:timeseries] = Dict()
    mod.ext[:variables] = Dict()
    mod.ext[:constraints] = Dict()
    mod.ext[:expressions] = Dict()

    # Sets
    mod.ext[:sets][:JH] = 1:data["General"]["nTimesteps"]
    mod.ext[:sets][:JZ] = 1:length(zones)
  
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

    for row in eachrow(participation_matrix)
        agent_name = row[:agent]
        participation_switch[agent_name] = Dict(z => row[Symbol(z)] for z in zones)
    end
    # data["General"]["participation_matrix"] = participation_switch
    mod.ext[:parameters][:participation_matrix] = participation_switch

    # Derating factors
    derating_switch = Dict{String, Dict{String, Float64}}()

    for row in eachrow(derating_factor)
        agent_name = row[:agent]
        derating_switch[agent_name] = Dict(z => row[Symbol(z)] for z in zones)
    end
    # data["General"]["derating_factor"] = derating_switch
    mod.ext[:parameters][:derating_factor] = derating_switch

    mod.ext[:parameters][:w] = ts[!, :weights][1:data["General"]["nTimesteps"]]
    
    # interconnectors
    mod.ext[:parameters][:TCONNECT] = [(Symbol(a),Symbol(b)) for (a,b) in data["Network"]["TCONNECT"]]
    # Node -> zone mapping
    zone_map = Dict{Symbol,Vector{Symbol}}()

    for (z,nodes) in  data["Network"]["ZoneMap"]
        zone_map[ Symbol(z) ] = Symbol.(nodes)
    end
    mod.ext[:parameters][:zone_map] = zone_map

    # list of nodes
    node_syms = reduce(vcat, values(zone_map))
    mod.ext[:parameters][:nodes] = node_syms

    # set of nodes
    mod.ext[:sets][:JN] = 1:length(node_syms)
    data["General"]["nNodes"] = length(node_syms)
    # node → zone lookup
    zone_of = Dict{Symbol,Symbol}()
    for (Z,N) in zone_map
        for n in N
            zone_of[n] = Z
        end
    end

    mod.ext[:parameters][:zone_of] = zone_of

    zone_syms = Symbol.(sort(zones))
    
    mod.ext[:parameters][:zone_syms] = zone_syms

    zone_idx = Dict(z=>i for (i,z) in enumerate(zone_syms))
    mod.ext[:parameters][:zone_of_idx] = [zone_idx[zone_of[n]] for n in node_syms]

    # Network topology
    mod.ext[:sets][:JL] = 1:nrow(lines)
    # Line IDs
    mod.ext[:parameters][:lines] = lines[!, :line_id]

    mod.ext[:parameters][:BRANCHES]  = [(Symbol(lines.from[i]), Symbol(lines.to[i]), lines.Fmax[i]) for i in 1:nrow(lines)]   # [(:n1, :n2, 4000), (:n2, :n4, 400), ...]

    mod.ext[:parameters][:zonal_PTDF] = Matrix(ptdf[!, zones])

    node_names = String.(mod.ext[:parameters][:nodes])
    mod.ext[:parameters][:nodal_PTDF] = Matrix(nodal_ptdf[:, node_names])

    mod.ext[:parameters][:G_nodal] =  zeros(data["General"]["nTimesteps"], length(node_syms))
    mod.ext[:parameters][:D_nodal] =  zeros(data["General"]["nTimesteps"], length(node_syms))
    mod.ext[:parameters][:Y_nodal] =  zeros(data["General"]["nTimesteps"], length(node_syms))


    # capacity market parameters
    mod.ext[:parameters][:coupling] = data["Network"]["coupling"]
    mod.ext[:parameters][:MEC] = data["Network"]["MEC"]
    

    return mod, agents
end