function solve_capacityIC_agent!(mod::Model)
    
    # Sets
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]           # scarcity scenarios
    
    # Parameters
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    # zone_syms = mod.ext[:parameters][:zone_syms]
    # BRANCHES = mod.ext[:parameters][:BRANCHES]
    # Fmax = [BRANCHES[jl][3] for jl in JL]
    # nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    # d_scarcity = mod.ext[:parameters][:d_scarcity]
    CapCM_zonal = mod.ext[:parameters][:CapCM_zonal]
    # Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal]
    rc = mod.ext[:parameters][:reserve_cost]
    coupling = mod.ext[:parameters][:coupling]

    y_bar_nodal = mod.ext[:parameters][:y_bar_nodal]
    
    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar]
    λ_CM = mod.ext[:parameters][:λ_CM]
    ρ_CM = mod.ext[:parameters][:ρ_CM]
    
    # Variables
    cap_cm = mod.ext[:variables][:cap_cm] 
    # r_cm = mod.ext[:variables][:r_cm]
    # flow_cm  = mod.ext[:variables][:flow_cm]
    s_cm = mod.ext[:variables][:s_cm]
    CapCM_nodal = mod.ext[:variables][:CapCM_nodal]
    g_scar = mod.ext[:variables][:g_scar]


    # Objective function
    mod.ext[:objective] = @objective(mod, Min,
        - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
        + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
        + sum(rc * s_cm[jn] for jn in JN)
        )

    # constrain cap_cm to zero
    # @constraint(mod, [jz in JZ], cap_cm[jz] == 0)


    if coupling == "FB"

        for js in JS, jn in JN
            delete(mod, mod.ext[:constraints][:cap_limit][js,jn])
        end
        mod.ext[:constraints][:cap_limit] =
            @constraint(mod, [js in JS, jn in JN], g_scar[js, jn] <= y_bar_nodal[jn] + s_cm[jn])

        for jz in JZ
            delete(mod, mod.ext[:constraints][:capacity_allocation][jz])
        end
        mod.ext[:constraints][:capacity_allocation] = @constraint(mod, [jz in JZ],
            sum(CapCM_nodal[jn] for jn in JN if zone_of_idx[jn] == jz) == CapCM_zonal[jz])

        for jn in JN
            delete(mod, mod.ext[:constraints][:capcm_limit][jn])
        end
        mod.ext[:constraints][:capcm_limit] = @constraint(mod, [jn in JN], CapCM_nodal[jn] <= y_bar_nodal[jn])

    elseif coupling == "ATC"
        
        TCONNECT = mod.ext[:parameters][:TCONNECT]
        # ATC = mod.ext[:parameters][:ATC]
        ex_cm = mod.ext[:variables][:ex_cm]

        for js in JS, jz in JZ
            delete(mod, mod.ext[:constraints][:cap_limit][js,jz])
        end
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [js in JS, jz in JZ], sum(g_scar[js, jn]  for jn in JN if zone_of_idx[jn] == jz) <= CapCM_zonal[jz] + sum(s_cm[jn] for jn in JN if zone_of_idx[jn] == jz))

        
        # ATC limits on exchanges
        # for js in JS, t in TCONNECT
        #     delete(mod, mod.ext[:constraints][:cap_cm_atc_limit][js,t])
        # end
        # mod.ext[:constraints][:cap_cm_atc_limit] = @constraint(mod, [js in JS, t in TCONNECT], 
        #     ATC[js][t][2] <= ex_cm[js,t] <= ATC[js][t][1])
    end

    optimize!(mod)

    return mod
end
