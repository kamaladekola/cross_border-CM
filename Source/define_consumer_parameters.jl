function define_consumer_parameters!(mod::Model, data::Dict,ts::DataFrame)
    # Parameters - note consumers are rescaled (total number of consumers x share of this type of consumer)
    mod.ext[:timeseries][:D] = data["totConsumers"]*data["Share"]*ts[!,data["D"]] # demand profile 
    mod.ext[:timeseries][:PV] = data["totConsumers"]*data["Share"]*data["PV_cap"]*ts[!,data["PV_AF"]] # pv profile 

    return mod
end