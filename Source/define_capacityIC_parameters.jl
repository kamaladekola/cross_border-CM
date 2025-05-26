function define_capacityIC_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, scarcity::DataFrame)

    mod.ext[:sets][:JS] = 1:nrow(scarcity)

    node_syms = mod.ext[:parameters][:nodes]

    node_shares = Dict{String,Dict{Symbol,Float64}}()
    for zone in zones
        cons_conf = get(data, zone, Dict{Any,Any}())
        shares = Dict{Symbol,Float64}()
        ns = get(cons_conf, "NodeShare", get(cons_conf, :NodeShare, Dict{Any,Any}()))
        for (n_str, s) in ns
            shares[Symbol(n_str)] = s
        end
        node_shares[zone] = shares
    end

    M = zeros(length(zones), length(node_syms))
    for (jz, zone) in enumerate(zones)
        for (jn, node) in enumerate(node_syms)
            M[jz, jn] = get( node_shares[zone], node, 0.0 )
        end
    end

    sc = copy(scarcity)
    rename!(sc, names(sc)[1] => :scenario)
    scarcity_zonal = select(sc, zones) |> Matrix
    mod.ext[:parameters][:d_scarcity] = scarcity_zonal * M

    # initialize the CM accumulators
    # mod.ext[:parameters][:CapCM_zonal] = zeros(length(zones))
    mod.ext[:parameters][:CapCM_nodal] = zeros(data["nNodes"])

    return mod
end