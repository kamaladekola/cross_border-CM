function define_generator_parameters!(mod::Model, data::Dict,ts::DataFrame, af::DataFrame)
    # Parameters 
    mod.ext[:parameters][:A] = data["a"]
    mod.ext[:parameters][:B] = data["b"]
    mod.ext[:parameters][:C] = data["C"]

    # investment parameters
    mod.ext[:parameters][:I] = data["I"]
    mod.ext[:parameters][:max_cap]  = data["max_cap"]

    mod.ext[:parameters][:σ_CM] = data["sigmaCM"]


   # Availability factors
    if haskey(data,"AF")
        mod.ext[:timeseries][:AF] = af[!,data["AF"]]
    else
        mod.ext[:timeseries][:AF] = ones(data["nTimesteps"])
    end 

    return mod
end