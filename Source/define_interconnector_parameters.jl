function define_interconnector_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame)
    
    # Store parameters in model
    mod.ext[:parameters][:lines] = ptdf[!, :line_id]
    mod.ext[:parameters][:RAM] = ptdf[!, :RAM]
    mod.ext[:parameters][:PTDF] =  Matrix(ptdf[!, zones])
    

    
    return mod
end