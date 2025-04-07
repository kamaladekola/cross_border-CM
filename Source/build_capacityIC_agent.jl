function build_capacityIC_agent!(mod::Model)
    JZ = mod.ext[:sets][:JZ]
    JL = mod.ext[:sets][:JL]

    RAM = mod.ext[:parameters][:RAM] 
    PTDF = mod.ext[:parameters][:PTDF]
    np_max = mod.ext[:parameters][:np_max]

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar] # ADMM consensus variables
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]

    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], base_name = "zonal capacity net position")

    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
)

    # Constraints
    mod.ext[:constraints][:CapacityNP] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0) # causes issue when only some zones are active in the CM
    mod.ext[:constraints][:CapacityIC] = @constraint(mod, [jl in JL], -RAM[jl] <= sum(PTDF[jl, jz] * cap_cm[jz] for jz in JZ ) <= RAM[jl])

    # one constraint for each scarcity event
    for (js, row) in enumerate(eachrow(np_max))
        mod.ext[:constraints][Symbol("MaxNP_$js")] = @constraint(mod, [jz in JZ], -row[jz] <= cap_cm[jz] <= row[jz])
    end
    return mod
end
