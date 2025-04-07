function build_interconnector_agent!(mod::Model)
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    JL = mod.ext[:sets][:JL]

    W = mod.ext[:parameters][:w]

    RAM = mod.ext[:parameters][:RAM] 
    PTDF = mod.ext[:parameters][:PTDF]


    λ_all = mod.ext[:parameters][:λ_all]  # λ_EOM[t,z]

    # ADMM consensus variables and penalty
    g_bar_all = mod.ext[:parameters][:g_bar_all]  # ḡ[t,z]
    ρ_all = mod.ext[:parameters][:ρ_all]

    # net-position variable g[t,z] --> positive => import
    g = mod.ext[:variables][:g] = @variable(mod, g[jh=JH,jz=JZ], base_name = "netposition")

    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    - sum(W[jh] * λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum(W[jh] * (ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 
    for jh in JH, jz in JZ))

    # Constraints
    # Net position constraint -->  sum_z g[t,z] = 0 for each time t
    mod.ext[:constraints][:netposition] = @constraint(mod, [jh=JH], sum(g[jh, jz] for jz in JZ) == 0)
    mod.ext[:constraints][:FBMC] = @constraint(mod,[jh=JH, jl=JL], -RAM[jl] <= sum(PTDF[jl, jz] * g[jh, jz] for jz in JZ) <= RAM[jl])
    
    return mod
end
