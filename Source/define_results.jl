function define_results!(data::Dict,results::Dict,ADMM::Dict,agents::Dict) 
    results["g"] = Dict()
    for m in agents[:eom]
        results["g"][m] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
        push!(results["g"][m],zeros(data["nTimesteps"]))
    end
 
    results["λ"] = Dict()
    results[ "λ"]["EOM"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results[ "λ"]["EOM"],zeros(data["nTimesteps"]))
   
    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["EOM"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["EOM"],zeros(data["nTimesteps"]))
  
    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = Dict()
    ADMM["Residuals"]["Primal"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["EOM"],0)

    ADMM["Residuals"]["Dual"] = Dict()
    ADMM["Residuals"]["Dual"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["EOM"],0)
    

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["EOM"] = data["epsilon"] 

    ADMM["ρ"] = Dict()
    ADMM["ρ"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["EOM"],data["rho_EOM"])

    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0
    
    return results, ADMM
end