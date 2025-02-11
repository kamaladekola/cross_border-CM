function define_consumer_parameters!(mod::Model, data::Dict, load::DataFrame)
    # Parameters - note consumers are rescaled (total number of consumers x share of this type of consumer)
    mod.ext[:timeseries][:D] = load[!, Symbol(data["D"])][1:data["nTimesteps"]] # demand profile 
    mod.ext[:parameters][:WTP] = data["WTP"] # value of lost load
    mod.ext[:parameters][:ela] = data["ela"] # fraction of demand that is elastic
    return mod
end