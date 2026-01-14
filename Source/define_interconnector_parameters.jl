function define_interconnector_parameters!(mod::Model, data::Dict, zones::Vector{String}, ptdf::DataFrame, nodal_ptdf::DataFrame, lines::DataFrame, weights)
    
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    
    W = JuMP.Containers.DenseAxisArray(
        [Float64(weights[jh, Symbol(zones[jz])]) for jh in JH, jz in JZ],
        JH, JZ
    )
    
    mod.ext[:parameters][:weight] = W

    @assert haskey(mod.ext[:sets], :JL) "JL not set; call define_common_parameters! first"
    @assert haskey(mod.ext[:parameters], :zonal_PTDF)  "zonal_PTDF missing"
    @assert haskey(mod.ext[:parameters], :nodal_PTDF)  "nodal_PTDF missing"
    @assert haskey(mod.ext[:parameters], :BRANCHES)    "BRANCHES missing"
    @assert haskey(mod.ext[:parameters], :nodes)       "nodes missing"

    return mod
end