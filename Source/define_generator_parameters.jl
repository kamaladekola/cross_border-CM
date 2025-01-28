function define_generator_parameters!(mod::Model, data::Dict,ts::DataFrame)
    # Parameters 
    mod.ext[:parameters][:A] = data["a"]
    mod.ext[:parameters][:B] = data["b"]
    mod.ext[:parameters][:C] = data["C"]

   # Availability factors
    if haskey(data,"AF")
        mod.ext[:timeseries][:AC] = data["C"]*ts[!,data["AF"]]  
    else
        mod.ext[:timeseries][:AC] = data["C"]*ones(data["nTimesteps"]) 
    end 
   
    return mod
end