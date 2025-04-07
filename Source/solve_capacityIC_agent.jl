function solve_capacityIC_agent!(mod::Model)
    
    JZ = mod.ext[:sets][:JZ]

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar]
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]

    cap_cm = mod.ext[:variables][:cap_cm]

    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
)

    optimize!(mod);

    return mod
end