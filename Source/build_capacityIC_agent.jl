function build_capacityIC_agent!(mod::Model)
    # Sets
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]

    # Parameters
    coupling  = mod.ext[:parameters][:coupling]
    nodes = mod.ext[:parameters][:nodes]
    zone_of = mod.ext[:parameters][:zone_of]
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    zone_syms = mod.ext[:parameters][:zone_syms]
    TCONNECT = mod.ext[:parameters][:TCONNECT]
    BRANCHES = mod.ext[:parameters][:BRANCHES]
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    d_scarcity = mod.ext[:parameters][:d_scarcity]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]

    ATC  = mod.ext[:parameters][:MEC]
    # MEC = get_ATC()

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar] # ADMM consensus variables
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]

    # Variables
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], base_name = "zonal_capacity_netpos")
    flow_cm = mod.ext[:variables][:flow_cm] = @variable(mod, [js=JS, jl=JL], base_name = "flow_cm")
    ex_cm = mod.ext[:variables][:ex_cm] = @variable(mod, [js=JS, t in TCONNECT], base_name = "ex_cm")
    g_scarcity = mod.ext[:variables][:g_scarcity] = @variable(mod, [js=JS,jn=JN], lower_bound=0, base_name = "dispatch")     # check dispatch feasibility
    ens_cm = mod.ext[:variables][:ens_cm] = @variable(mod, [js=JS, jn=JN], lower_bound=0, base_name = "ens_cm")
    g_red_cm = mod.ext[:variables][:g_red_cm] = @variable(mod, g_red[js in JS, jn in JN], base_name = "g_red_cm")

    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
    + sum(((ens_cm[js,jn])^2 * 1e10) for js in JS, jn in JN)
    )


    # mod.ext[:constraints][:CapacityNP] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0) # capacity netposition

    # # scenario based dispatch to simulate different scarcity moments across zones
    # # mod.ext[:constraints][:dispatch] = @constraint(mod, [js in JS, jn in JN], g_scarcity[js, jn] == d_scarcity[js, jn] ) # include nodal clearing -> g_scarcity = D_scarcity (to ensure feasible dispatch)
    
    # mod.ext[:constraints][:dispatch_feasibility] = @constraint(mod, [js in JS, jn in JN], g_scarcity[js, jn] ≤ CapCM_nodal[jn])  # dispatched electricity across interconnector must be less than capacity offers.
    
    # # ############ Including exact projection constraints
    # # exchanges and flows
    # for (jl,(i,j,_)) in enumerate(BRANCHES), js in JS
    #     if (zone_of[i],zone_of[j]) in TCONNECT
    #         mod.ext[:constraints][Symbol("ex_cm$(js)_$(jl)")] = @constraint(mod, flow_cm[js,jl] == ex_cm[js, (zone_of[i],zone_of[j])])
    #     elseif (zone_of[j],zone_of[i]) in TCONNECT
    #         mod.ext[:constraints][Symbol("ex_cm$(js)_$(jl)")] = @constraint(mod, flow_cm[js,jl] == -ex_cm[js, (zone_of[j],zone_of[i])])
    #     end
    # end

    # # netposition = exports + imports (negative netposition --> export) 
    # mod.ext[:constraints][:cap_cm_netposition] = @constraint(mod, [js in JS, jz in JZ],
    #     cap_cm[jz] ==
    #       sum(ex_cm[js,t] for t in TCONNECT if t[2] == zone_syms[jz])
    #       - sum(ex_cm[js,t] for t in TCONNECT if t[1] == zone_syms[jz])
    #       )

    # if coupling == "FB"
    #     # nodal balance constraint: flow = generation + redispatch - demand --> for all nodes
    #     mod.ext[:constraints][:cap_cm_nodal_balance] = @constraint(mod, [js in JS, jl in JL],
    #         flow_cm[js,jl] == sum(nodal_PTDF[jl, jn] * (g_scarcity[js,jn] - d_scarcity[js,jn] + ens_cm[js,jn]) for jn in JN))

    #     # thermal limit
    #     mod.ext[:constraints][:thermal_limit] = @constraint(mod, [js in JS, jl in JL], -BRANCHES[jl][3] ≤ flow_cm[js,jl] ≤ BRANCHES[jl][3])

    #     ################# debug ###########################
    #     # mod.ext[:constraints][:redispatch_limit] = @constraint(mod, [js in JS, jn in JN], 
    #     #     0 <= g_scarcity[js,jn] + g_red_cm[js,jn] <= CapCM_nodal[jn]) 

    #     # mod.ext[:constraints][:redispatch_balance] = @constraint(mod, [js in JS, jz in JZ],
    #     #       sum(g_red_cm[js,jn] for jn in JN if zone_of_idx[jn] == jz) == 0)

    # elseif coupling == "ATC"     # each border is independently constrained by ATC 
    #    ATC = mod.ext[:parameters][:ATC]    # Dict(A,B) 
    #    mod.ext[:constraints][:cap_cm_atc_limit] =     @constraint(mod, [s in JS, t in TCONNECT],
    #     -ATC[t] ≤ ex_cm[s,t] ≤ ATC[t])
    # end

    return mod
end

# dscarcity as ratios of total capacity offered
