function define_generator_parameters!(mod::Model, data::Dict, weights::DataFrame, af::DataFrame, zone::SubString{String})
    # Parameters 
    mod.ext[:parameters][:A] = data["a"]
    mod.ext[:parameters][:B] = data["b"]
    mod.ext[:parameters][:C] = data["C"]

    mod.ext[:parameters][:w] = weights[!, Symbol(zone)][1:data["nTimesteps"]]
    # investment parameters
    mod.ext[:parameters][:I] = data["I"]
    mod.ext[:parameters][:max_cap]  = data["max_cap"]

    mod.ext[:parameters][:σ_CM] = data["sigmaCM"]

    # Nodal allocation
    nodes = mod.ext[:parameters][:nodes]  # Vector{Symbol}
    raw_ns = Dict{Symbol,Float64}()
    for (node_str, share) in data["NodeShare"]
        raw_ns[ Symbol(node_str) ] = float(share)
    end
    node_share_vec = [ get(raw_ns, n, 0.0) for n in nodes ]  # length |N|
    mod.ext[:parameters][:node_share_vec] = node_share_vec

   # Availability factors
    if haskey(data,"AF")
        mod.ext[:timeseries][:AF] = af[!,data["AF"]]
    else
        mod.ext[:timeseries][:AF] = ones(data["nTimesteps"])
    end 
    

    
    return mod
end