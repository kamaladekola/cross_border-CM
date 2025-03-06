function solve_interconnector_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]



    λ_all = mod.ext[:parameters][:λ_all]  # λ_EOM[t,z]

    # ADMM consensus variables and penalty
    g_bar_all = mod.ext[:parameters][:g_bar_all]  # ḡ[t,z]
    ρ_all = mod.ext[:parameters][:ρ_all]
    # net position: positive => import

    g = mod.ext[:variables][:g]

    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    -sum(λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum((ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 
    for jh in JH, jz in JZ))



    optimize!(mod);

    return mod
end
