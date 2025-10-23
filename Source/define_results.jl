function define_results!(data::Dict, results::Dict, ADMM::Dict, agents::Dict, zones::Vector)
    # Store generation results per agent
    results["g"] = Dict()
    for m in agents[:eom]
        if m == "NetworkManager"
            results["g"][m] = CircularBuffer{Matrix{Float64}}(data["CircularBufferSize"])
            push!(results["g"][m], zeros(data["nTimesteps"], length(zones)))
        else
            results["g"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
            push!(results["g"][m], zeros(data["nTimesteps"]))
        end
    end

    # Store installed capacity results for generators
    results["y"] = Dict()
    results["y_nodal"] = Dict()
    results["CapCM_nodal"] = Dict()

    for m in agents[:Gen]
        results["y"][m] = CircularBuffer{Float64}(data["CircularBufferSize"])
        results["y_nodal"][m] = CircularBuffer{Matrix{Float64}}(data["CircularBufferSize"])
        results["CapCM_nodal"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
        push!(results["y"][m], 0.0)
        push!(results["y_nodal"][m], zeros(data["nTimesteps"], data["nNodes"]))
        push!(results["CapCM_nodal"][m], zeros(data["nNodes"]))
    end

    # Store consumer results
    results["Cons"] = Dict()
    results["Cons"]["inelastic_demand"] = Dict()
    results["Cons"]["elastic_demand"] = Dict()
    results["Cons"]["ENS"] = Dict()
    results["Cons"]["d_nodal"] = Dict()
    for m in agents[:Cons]
       results["Cons"]["inelastic_demand"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
       results["Cons"]["elastic_demand"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
       results["Cons"]["ENS"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
       results["Cons"]["d_nodal"][m] = CircularBuffer{Matrix{Float64}}(data["CircularBufferSize"])

       push!(results["Cons"]["inelastic_demand"][m], zeros(data["nTimesteps"]))
       push!(results["Cons"]["elastic_demand"][m], zeros(data["nTimesteps"]))
       push!(results["Cons"]["ENS"][m], zeros(data["nTimesteps"]))
       push!(results["Cons"]["d_nodal"][m], zeros(data["nTimesteps"], data["nNodes"]))
    end

    # Store capacity offered by agents in the capacity market
    results["cap_cm"] = Dict()
    for m in agents[:cm]
        results["cap_cm"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
        push!(results["cap_cm"][m], zeros(length(zones)))
    end

    # Store y_bar_nodal results
    results["y_bar_nodal"] = Dict()
    for m in agents[:IC]
        results["y_bar_nodal"][m] = CircularBuffer{Vector{Float64}}(data["CircularBufferSize"])
        push!(results["y_bar_nodal"][m], zeros(data["nNodes"]))
    end

    # Prices for each zone in the EOM market
    results["λ"] = Dict()
    results["λ"]["EOM"] = Dict(z => CircularBuffer{Vector{Float64}}(data["CircularBufferSize"]) for z in zones)
    results["λ"]["CM"] = Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
    for z in zones
        push!(results["λ"]["EOM"][z], zeros(data["nTimesteps"]))
        push!(results["λ"]["CM"][z], 0.0)
    end

    # Imbalances per zone
    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["EOM"] = Dict(z => CircularBuffer{Vector{Float64}}(data["CircularBufferSize"]) for z in zones)
    ADMM["Imbalances"]["CM"] = Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
    for z in zones
        push!(ADMM["Imbalances"]["EOM"][z], zeros(data["nTimesteps"]))
        push!(ADMM["Imbalances"]["CM"][z], 0.0)
    end

    # Residuals per zone (Primal and Dual)
    ADMM["Residuals"] = Dict(
        "Primal" => Dict(
            "EOM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones),
            "CM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
        ),
        "Dual" => Dict(
            "EOM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones),
            "CM" => Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
        )
    )
    for z in zones
        push!(ADMM["Residuals"]["Primal"]["EOM"][z], 0.0)
        push!(ADMM["Residuals"]["Dual"]["EOM"][z], 0.0)
        push!(ADMM["Residuals"]["Primal"]["CM"][z], 0.0)
        push!(ADMM["Residuals"]["Dual"]["CM"][z], 0.0)
    end
    
    # Tolerance for EOM and CM
    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["EOM"] = data["epsilon"]
    ADMM["Tolerance"]["CM"] = data["epsilon"]

    # Initialize per-zone ρ (rho) values
    ADMM["ρ"] = Dict()
    ADMM["ρ"]["EOM"] = Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
    ADMM["ρ"]["CM"] = Dict(z => CircularBuffer{Float64}(data["CircularBufferSize"]) for z in zones)
    for z in zones
        push!(ADMM["ρ"]["EOM"][z], data["rho_EOM"])
        push!(ADMM["ρ"]["CM"][z], data["rho_CM"])
    end

    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0
    
    return results, ADMM
end
