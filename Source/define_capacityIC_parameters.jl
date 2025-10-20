function define_capacityIC_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, scarcity::DataFrame)
    # Define scarcity scenarios set
    mod.ext[:sets][:JS] = 1:nrow(scarcity)
    JS = mod.ext[:sets][:JS]
    
    # Get node symbols from common parameters
    node_syms = mod.ext[:parameters][:nodes]
    # Build nodal distribution matrix for consumers
    node_shares = Dict{String, Dict{Symbol,Float64}}()
    for zone in zones
        cons_conf = data["Consumers"][zone]
        raw_ns = get(cons_conf, "NodeShare", Dict())
        node_shares[zone] = Dict(Symbol(k) => v for (k,v) in raw_ns)
    end
    
    # Node share matrix: zones × nodes
    M = [get(node_shares[zone], node, 0.0) for zone in zones, node in node_syms]
    

    # sc = copy(scarcity)
    # rename!(sc, names(sc)[1] => :scenario)
    # scarcity_zonal = select(sc, zones) |> Matrix
    # mod.ext[:parameters][:d_scarcity] = scarcity_zonal * M
    
    scarcity_matrix = Matrix{Float64}(scarcity[:, zones])
    mod.ext[:parameters][:scarcity_matrix] = scarcity_matrix
    
    # Set capacity market parameters
    mod.ext[:parameters][:CapCM_zonal] = zeros(length(zones))
    capacity_demand_zonal = [get(data["CM"][zone], "capacity_target", 0.0) for zone in zones]
    mod.ext[:parameters][:Cap_Demand_nodal] = vec(M' * capacity_demand_zonal)
    
    
    # Store TCONNECT from network configuration
    if haskey(data["Network"], "TCONNECT")
        # Convert string connections to symbol tuples
        mod.ext[:parameters][:TCONNECT] = [(Symbol(t[1]), Symbol(t[2])) for t in data["Network"]["TCONNECT"]]
    else
        # Create default connections between all zones if not specified
        zone_syms = Symbol.(zones)
        mod.ext[:parameters][:TCONNECT] = [(z1, z2) for z1 in zone_syms for z2 in zone_syms if z1 != z2]
    end
    
    # ATC parameter:
    mod.ext[:parameters][:ATC] = Dict{Tuple{Symbol,Symbol}, Tuple{Float64,Float64}}()

    # for t in mod.ext[:parameters][:TCONNECT]
    #     mod.ext[:parameters][:ATC][t] = (3000.0, -3000.0)
    # end

    if (:A, :C) in mod.ext[:parameters][:TCONNECT]
        mod.ext[:parameters][:ATC][(:A, :C)] = (2000.0, -2000.0)
    end

    if (:C, :B) in mod.ext[:parameters][:TCONNECT]
        mod.ext[:parameters][:ATC][(:C, :B)] = (2000.0, -2000.0)
    end

    if (:B, :A) in mod.ext[:parameters][:TCONNECT]
        mod.ext[:parameters][:ATC][(:B, :A)] = (1000.0, -1000.0)
    end

    return mod
end