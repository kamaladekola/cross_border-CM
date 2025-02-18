function define_consumer_parameters!(mod::Model, data::Dict, load::DataFrame)

    mod.ext[:timeseries][:D] = load[!, Symbol(data["D"])][1:data["nTimesteps"]] # demand profile 
    mod.ext[:parameters][:WTP] = data["WTP"] # value of lost load
    mod.ext[:parameters][:ela] = data["ela"] # fraction of demand that is elastic
    mod.ext[:parameters][:D_max] = maximum(load[!, Symbol(data["D"])][1:data["nTimesteps"]])
    mod.ext[:parameters][:σ_CM] = data["sigmaCM"]

    return mod
end