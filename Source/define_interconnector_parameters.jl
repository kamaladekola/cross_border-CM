function define_interconnector_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, nodal_ptdf::DataFrame, lines::DataFrame)
    
    mod.ext[:sets][:JL] = 1:nrow(lines)
    # Line IDs
    mod.ext[:parameters][:lines] = lines[!, :line_id]

    BRANCHES = [(Symbol(lines.from[i]), Symbol(lines.to[i]), lines.Fmax[i]) for i in 1:nrow(lines)]   # [(:n1, :n2, 4000), (:n2, :n4, 400), ...]

    mod.ext[:parameters][:BRANCHES] = BRANCHES

    mod.ext[:parameters][:zonal_PTDF] = Matrix(ptdf[!, zones])

    node_names = String.(mod.ext[:parameters][:nodes])
    mod.ext[:parameters][:nodal_PTDF] = Matrix(nodal_ptdf[:, node_names])


    node_syms = mod.ext[:parameters][:nodes]

    mod.ext[:parameters][:G_nodal] =  zeros(data["nTimesteps"], length(node_syms))
    mod.ext[:parameters][:D_nodal] =  zeros(data["nTimesteps"], length(node_syms))
    mod.ext[:parameters][:Y_nodal] =  zeros(data["nTimesteps"], length(node_syms))

    return mod
end

