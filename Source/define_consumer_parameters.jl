function define_consumer_parameters!(mod::Model, data::Dict, load::DataFrame, weights::DataFrame)

    mod.ext[:parameters][:w] = weights[!, Symbol(data["D"])][1:data["nTimesteps"]]
    mod.ext[:timeseries][:D] = load[!, Symbol(data["D"])][1:data["nTimesteps"]] # demand profile 
    mod.ext[:parameters][:WTP] = data["WTP"] # value of lost load
    mod.ext[:parameters][:ela] = data["ela"] # fraction of demand that is elastic
    mod.ext[:parameters][:D_max] = maximum(load[!, Symbol(data["D"])][1:data["nTimesteps"]])
    mod.ext[:parameters][:σ_CM] = data["sigmaCM"]

    # CM parameters
    mod.ext[:parameters][:CD] = data["capacity_target"]          # Get capacity target for zone Z
    mod.ext[:parameters][:WTP_CM] = data["price_target"]         # Willingness to pay for capacity in the CM
    mod.ext[:parameters][:CD_margin] = data["capacity_margin"]   # Minimum willingness to pay for capacity in the CM

    # Nodal allocation
    nodes = mod.ext[:parameters][:nodes]  # Vector{Symbol}
    raw_ns = Dict{Symbol,Float64}()
    for (node_str, share) in data["NodeShare"]
        raw_ns[ Symbol(node_str) ] = float(share)
    end

    node_share_vec = [ get(raw_ns, n, 0.0) for n in nodes ]  # length |N|
    mod.ext[:parameters][:node_share_vec] = node_share_vec

    return mod
end
