function define_capacityIC_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, RAM_scen::DataFrame)

    # Store parameters in model
    
    # mod.ext[:parameters][:lines] = ptdf[:, :line_id]
    line_ids = ptdf[:, :line_id]
    mod.ext[:parameters][:RAM] = ptdf[!, :RAM]
    mod.ext[:parameters][:PTDF] =  Matrix(ptdf[!, zones])
    # mod.ext[:parameters][:RAM_scenarios] = hcat([RAM_scen[:, Symbol(l)] for l in line_ids]...) |> Matrix
    mod.ext[:parameters][:coupling] = data["coupling"]
    mod.ext[:parameters][:MEC] = data["MEC"]
    mod.ext[:parameters][:from_zone] = data["from_zone"]
    mod.ext[:parameters][:to_zone] = data["to_zone"]
    # mod.ext[:sets][:border_names] = ["AB", "BC", "CA"]
    
    return mod
end