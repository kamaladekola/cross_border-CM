function solve_interconnector_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    JL = mod.ext[:sets][:JL]
    JN = mod.ext[:sets][:JN]
    W = mod.ext[:parameters][:w]

    G_nodal = mod.ext[:parameters][:G_nodal]         # |H|×|N|
    D_nodal = mod.ext[:parameters][:D_nodal]         # |H|×|N|
    Y_nodal = mod.ext[:parameters][:Y_nodal]         # |H|×|N|
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]      # |L|×|N| nodal PTDF

    λ_all = mod.ext[:parameters][:λ_all]  # λ_EOM[t,z]

    # ADMM consensus variables and penalty
    g_bar_all = mod.ext[:parameters][:g_bar_all]  # ḡ[t,z]
    ρ_all = mod.ext[:parameters][:ρ_all]
    # net position: positive => import

    g = mod.ext[:variables][:g]
    g_red = mod.ext[:variables][:redispatch]
    flow = mod.ext[:variables][:flow]
    ens_IC = mod.ext[:variables][:ens_IC]

    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    - sum(W[jh] * λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum(W[jh] * (ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 for jh in JH, jz in JZ)
    + sum(W[jh] * ((ens_IC[jh,jn])^2 * 1e5) for jh in JH, jn in JN)
    )


    for jh in JH, jl in JL
        delete(mod, mod.ext[:constraints][:nodal_balance][jh,jl])
    end
    mod.ext[:constraints][:nodal_balance] = @constraint(mod, [jh in JH, jl in JL],
        flow[jh,jl] == sum(nodal_PTDF[jl, jn] * (G_nodal[jh,jn] + g_red[jh,jn] - (D_nodal[jh,jn] + ens_IC[jh,jn])) for jn in JN))


        
    for jh in JH, jn in JN
        delete(mod, mod.ext[:constraints][:redispatch_limit][jh,jn])
    end
    mod.ext[:constraints][:redispatch_limit] = @constraint(mod, [jh in JH, jn in JN], 
    0 <= G_nodal[jh,jn] + g_red[jh,jn] + ens_IC[jh,jn] <= Y_nodal[jh,jn])

    optimize!(mod);

    return mod
end


