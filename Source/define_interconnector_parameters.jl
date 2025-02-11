function define_interconnector_parameters!(mod::Model, data::Dict, ts::DataFrame)
    IC = data["f_cap"]
    IC_r = data["r_cap"]
    
    mod.ext[:parameters][:IC] = IC
    mod.ext[:parameters][:IC_r] = IC_r
    mod.ext[:parameters][:λ1] = zeros(data["nTimesteps"])
    mod.ext[:parameters][:λ2] = zeros(data["nTimesteps"])
    
    return mod
end