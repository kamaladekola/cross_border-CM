function define_capacityIC_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, scarcity::DataFrame)

    mod.ext[:sets][:JS] = 1:nrow(scarcity)
    JS = mod.ext[:sets][:JS]

    node_syms = mod.ext[:parameters][:nodes]

    node_shares = Dict{String, Dict{Symbol,Float64}}()
    for zone in zones
        # data["Consumers"][zone] still contains the NodeShare block
        cons_conf = data["Consumers"][zone]
        raw_ns = get(cons_conf, "NodeShare", Dict())  # Dict("n1"=>0.51, …)
        node_shares[zone] = Dict(Symbol(k) => v for (k,v) in raw_ns)
    end


    M = [get(node_shares[zone], node, 0.0) for zone in zones, node in node_syms]

    sc = copy(scarcity)
    rename!(sc, names(sc)[1] => :scenario)
    scarcity_zonal = select(sc, zones) |> Matrix
    mod.ext[:parameters][:d_scarcity] = scarcity_zonal * M

    # initialize the CM accumulators
    # mod.ext[:parameters][:CapCM_zonal] = zeros(length(zones))
    TCONNECT = mod.ext[:parameters][:TCONNECT]
    mod.ext[:parameters][:CapCM_nodal] = zeros(data["General"]["nNodes"])
    mod.ext[:parameters][:ATC] = Dict(js => Dict(t => (0.0, 0.0) for t in TCONNECT)  for js in JS)

    capacity_demand_zonal = [get(data["CM"][zone], "capacity_target", 0.0) for zone in zones]
    
    mod.ext[:parameters][:Cap_Demand_nodal] = vec(M' * capacity_demand_zonal)
    # mod.ext[:parameters][:Cap_Demand_nodal] = zeros(data["General"]["nNodes"]) # debugging

    return mod
end