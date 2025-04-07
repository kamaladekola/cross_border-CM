function define_capacityIC_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, np_max::DataFrame)

    # Store parameters in model
    mod.ext[:parameters][:lines] = ptdf[!, :line_id]
    mod.ext[:parameters][:RAM] = ptdf[!, :RAM]
    mod.ext[:parameters][:PTDF] =  Matrix(ptdf[!, zones])
    mod.ext[:parameters][:np_max] = Matrix(np_max[!, zones])
    
    return mod
end