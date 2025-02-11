function define_results!(data::Dict, results::Dict, ADMM::Dict, agents::Dict, zones::Vector)
    # Store generation results per agent
    results["g"] = Dict()
    for m in agents[:eom]
        results["g"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
        push!(results["g"][m], zeros(data["nTimesteps"]))
    end

    # Prices for each zone in the EOM market
    results["λ"] = Dict()
    results["λ"]["EOM"] = Dict(z => CircularBuffer{Vector{Float64}}(data["CircularBufferSize"]) for z in zones)
    for z in zones
        push!(results["λ"]["EOM"][z], zeros(data["nTimesteps"]))
    end

    # Imbalances per zone
    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["EOM"] = Dict(z => CircularBuffer{Vector{Float64}}(data["CircularBufferSize"]) for z in zones)
    for z in zones
        push!(ADMM["Imbalances"]["EOM"][z], zeros(data["nTimesteps"]))
    end

    # Residuals per zone (Primal and Dual)
    ADMM["Residuals"] = Dict(
        "Primal" => Dict("EOM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)),
        "Dual"   => Dict("EOM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones))
    )
    for z in zones
        push!(ADMM["Residuals"]["Primal"]["EOM"][z], 0.0)
        push!(ADMM["Residuals"]["Dual"]["EOM"][z], 0.0)
    end
    
    # Tolerance for EOM (assumed same for all zones)
    ADMM["Tolerance"] = Dict("EOM" => data["epsilon"])

    # Initialize per-zone ρ (rho) values
    ADMM["ρ"] = Dict("EOM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones))
    for z in zones
        push!(ADMM["ρ"]["EOM"][z], data["rho_EOM"])
    end

    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0
    
    return results, ADMM
end
