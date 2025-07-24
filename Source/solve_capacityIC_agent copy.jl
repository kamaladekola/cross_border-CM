function solve_capacityIC_agent!(mod::Model)
    
    JZ = mod.ext[:sets][:JZ]
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]
    JL = mod.ext[:sets][:JL]           # lines

    coupling    = mod.ext[:parameters][:coupling]
    d_scarcity  = mod.ext[:parameters][:d_scarcity]
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    cap_bar     = mod.ext[:parameters][:cap_bar]
    λ_CM        = mod.ext[:parameters][:λ_CM]
    ρ_CM        = mod.ext[:parameters][:ρ_CM]
    nodal_PTDF  = mod.ext[:parameters][:nodal_PTDF]
    Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal] # capacity demand per zone

    # cap_cm      = mod.ext[:variables][:cap_cm]
    cap_import = mod.ext[:variables][:cap_import]
    cap_export = mod.ext[:variables][:cap_export]
    cap_cm      = mod.ext[:expressions][:cap_cm] 
    flow_cm     = mod.ext[:variables][:flow_cm]
    ex_cm       = mod.ext[:variables][:ex_cm]
    g_scarcity  = mod.ext[:variables][:g_scarcity]
    # ens_cm      = mod.ext[:variables][:ens_cm]
    ens_pos     = mod.ext[:variables][:ens_pos]
    ens_neg     = mod.ext[:variables][:ens_neg]
    ens_cm      = mod.ext[:expressions][:ens_cm]
    # ens_cm      = mod.ext[:variables][:ens_cm]

    TCONNECT = mod.ext[:parameters][:TCONNECT]
    ATC = mod.ext[:parameters][:ATC]
    
    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #         d_scarcity[js,jn] * sum(CapCM_nodal[n] for n in JN))

    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #         d_scarcity[js, jn] * Cap_Demand_nodal[jn])

    # # Objective
    mod.ext[:objective] = @objective(mod, Min,
        - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ)
        + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
        + sum((ens_pos[js,jn] + ens_neg[js,jn]) * 1e10 for js in JS, jn in JN)
    )

    # mod.ext[:objective] = @objective(mod, Min, 0)

    # nodal capacity limits
    for js in JS, jn in JN
        delete(mod, mod.ext[:constraints][:capacity_limit][js, jn])
    end
    mod.ext[:constraints][:capacity_limit] = @constraint(mod, [js in JS, jn in JN],
        g_scarcity[js, jn] <= CapCM_nodal[jn]
    )

    # # zonal balance
    # for js in JS, jz in JZ
    #     delete(mod, mod.ext[:constraints][:zonal_balance][js, jz])
    # end
    # mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
    #     sum(g_scarcity[js, jn] for jn in JN if zone_of_idx[jn] == jz)
    #     + sum(ens_cm[js, jn] for jn in JN if zone_of_idx[jn] == jz) + cap_cm[jz]
    #     == sum(demand[js, jn] for jn in JN if zone_of_idx[jn] == jz)
    # )


    if coupling == "FB"
        # for js in JS, jl in JL
        #     delete(mod, mod.ext[:constraints][:cap_cm_nodal_balance][js, jl])
        # end
        # mod.ext[:constraints][:cap_cm_nodal_balance] = @constraint(mod, [js in JS, jl in JL],
        #     flow_cm[js, jl] == sum(nodal_PTDF[jl, jn] * (g_scarcity[js, jn] - demand[js, jn] + ens_cm[js, jn]) for jn in JN)
        # )

        # for js in JS
        #     delete(mod, mod.ext[:constraints][:cap_cm_global_balance][js])
        # end
        # mod.ext[:constraints][:cap_cm_global_balance] = @constraint(mod, [js in JS],
        #     0 == sum((g_scarcity[js,jn] - demand[js,jn] + ens_cm[js,jn]) for jn in JN))


    elseif coupling == "ATC"
        
        # for js in JS, jz in JZ
        #     delete(mod, mod.ext[:constraints][:zonal_balance][js,jz])
        # end
        # mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
        # sum(g_scarcity[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + sum(ens_cm[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + cap_cm[jz] == sum(demand[js,jn] for jn in JN if zone_of_idx[jn] == jz))

        # for js in JS, t in TCONNECT
        #     delete(mod, mod.ext[:constraints][:cap_cm_atc_limit][js,t])
        # end
        mod.ext[:constraints][:cap_cm_atc_limit] = @constraint(mod, [js in JS, t in TCONNECT], 
            ATC[js][t][2] <= ex_cm[js,t] <= ATC[js][t][1])
         

    end

    optimize!(mod)

    return mod
end

# -------------------------------------------------------------------------------------
# old version
# -------------------------------------------------------------------------------------


function solve_capacityIC_agent!(mod::Model)
    
    JZ = mod.ext[:sets][:JZ]
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]
    JL = mod.ext[:sets][:JL]           # lines

    coupling    = mod.ext[:parameters][:coupling]
    d_scarcity  = mod.ext[:parameters][:d_scarcity]
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    cap_bar     = mod.ext[:parameters][:cap_bar]
    λ_CM        = mod.ext[:parameters][:λ_CM]
    ρ_CM        = mod.ext[:parameters][:ρ_CM]
    nodal_PTDF  = mod.ext[:parameters][:nodal_PTDF]
    Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal] # capacity demand per zone

    cap_cm      = mod.ext[:variables][:cap_cm]
    flow_cm     = mod.ext[:variables][:flow_cm]
    ex_cm       = mod.ext[:variables][:ex_cm]
    g_scarcity  = mod.ext[:variables][:g_scarcity]
    # ens_cm      = mod.ext[:variables][:ens_cm]
    ens_pos     = mod.ext[:variables][:ens_pos]
    ens_neg     = mod.ext[:variables][:ens_neg]
    ens_cm      = mod.ext[:expressions][:ens_cm]
    # ens_cm      = mod.ext[:variables][:ens_cm]

    TCONNECT = mod.ext[:parameters][:TCONNECT]
    ATC = mod.ext[:parameters][:ATC]
    
    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #         d_scarcity[js,jn] * sum(CapCM_nodal[n] for n in JN))

    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #         d_scarcity[js, jn] * Cap_Demand_nodal[jn])

    # # Objective
    mod.ext[:objective] = @objective(mod, Min,
        - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ)
        + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
        + sum((ens_pos[js,jn] + ens_neg[js,jn]) * 1e10 for js in JS, jn in JN)
    )

    # mod.ext[:objective] = @objective(mod, Min, 0)

    # nodal capacity limits
    for js in JS, jn in JN
        delete(mod, mod.ext[:constraints][:capacity_limit][js, jn])
    end
    mod.ext[:constraints][:capacity_limit] = @constraint(mod, [js in JS, jn in JN],
        g_scarcity[js, jn] <= CapCM_nodal[jn]
    )


    # # zonal balance
    # for js in JS, jz in JZ
    #     delete(mod, mod.ext[:constraints][:zonal_balance][js, jz])
    # end
    # mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
    #     sum(g_scarcity[js, jn] for jn in JN if zone_of_idx[jn] == jz)
    #     + sum(ens_cm[js, jn] for jn in JN if zone_of_idx[jn] == jz) + cap_cm[jz]
    #     == sum(demand[js, jn] for jn in JN if zone_of_idx[jn] == jz)
    # )


    if coupling == "FB"
        # for js in JS, jl in JL
        #     delete(mod, mod.ext[:constraints][:cap_cm_nodal_balance][js, jl])
        # end
        # mod.ext[:constraints][:cap_cm_nodal_balance] = @constraint(mod, [js in JS, jl in JL],
        #     flow_cm[js, jl] == sum(nodal_PTDF[jl, jn] * (g_scarcity[js, jn] - demand[js, jn] + ens_cm[js, jn]) for jn in JN)
        # )

        # for js in JS
        #     delete(mod, mod.ext[:constraints][:cap_cm_global_balance][js])
        # end
        # mod.ext[:constraints][:cap_cm_global_balance] = @constraint(mod, [js in JS],
        #     0 == sum((g_scarcity[js,jn] - demand[js,jn] + ens_cm[js,jn]) for jn in JN))


    elseif coupling == "ATC"
        
        # for js in JS, jz in JZ
        #     delete(mod, mod.ext[:constraints][:zonal_balance][js,jz])
        # end
        # mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
        # sum(g_scarcity[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + sum(ens_cm[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + cap_cm[jz] == sum(demand[js,jn] for jn in JN if zone_of_idx[jn] == jz))

        # for js in JS, t in TCONNECT
        #     delete(mod, mod.ext[:constraints][:cap_cm_atc_limit][js,t])
        # end
        mod.ext[:constraints][:cap_cm_atc_limit] = @constraint(mod, [js in JS, t in TCONNECT], 
            ATC[js][t][2] <= ex_cm[js,t] <= ATC[js][t][1])
            # -1000 <= ex_cm[js,t] <= 1000) # debugging

    end

    optimize!(mod)

    return mod
end
