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
    Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal] # capacity demand per zone
    ATC = mod.ext[:parameters][:ATC]

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar] # ADMM consensus variables
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]

    # Variables
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], base_name = "zonal_capacity_netpos")
    flow_cm = mod.ext[:variables][:flow_cm] = @variable(mod, [js=JS, jl=JL], base_name = "flow_cm")
    ex_cm = mod.ext[:variables][:ex_cm] = @variable(mod, [js=JS, t in TCONNECT], base_name = "ex_cm")
    g_scarcity = mod.ext[:variables][:g_scarcity] = @variable(mod, [js=JS,jn=JN], lower_bound=0, base_name = "dispatch")     # check dispatch feasibility

    ens_pos = mod.ext[:variables][:ens_pos] = @variable(mod, [js=JS, jn=JN], lower_bound = 0)          # extra gen   (MW)
    ens_neg = mod.ext[:variables][:ens_neg] = @variable(mod, [js=JS, jn=JN], lower_bound = 0)           # shed load   (MW)

    ens_cm = mod.ext[:expressions][:ens_cm] = @expression(mod, [js=JS, jn=JN], ens_pos[js,jn] - ens_neg[js,jn])
    # ens_cm = mod.ext[:variables][:ens_cm] = @variable(mod, [js=JS, jn=JN], lower_bound=0, base_name = "ens_cm")

    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
    # + sum(((ens_cm[js,jn])^2 * 1e8) for js in JS, jn in JN)
    )

    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #         d_scarcity[js, jn] * sum(CapCM_nodal[n] for n in JN))

    demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
            d_scarcity[js, jn] * Cap_Demand_nodal[jn])

    mod.ext[:constraints][:CapacityNP] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0) # capacity netposition

    # dispatch scenarios
    mod.ext[:constraints][:capacity_limit] = @constraint(mod, [js in JS, jn in JN], g_scarcity[js, jn] <= CapCM_nodal[jn])  # dispatched electricity across interconnector must be less than capacity offers.
    



          
    # ## netposition = exports + imports (negative netposition --> export) 
    # mod.ext[:constraints][:cap_cm_netposition] = @constraint(mod, [js in JS, jz in JZ],
    #     cap_cm[jz] == 0) # set this for implict CM and no cm
    # netposition = exports + imports (negative netposition --> export) 



    if coupling == "FB"
        # nodal balance constraint: flow = generation + redispatch - demand --> for all nodes
        mod.ext[:constraints][:cap_cm_nodal_balance] = @constraint(mod, [js in JS, jl in JL],
            flow_cm[js,jl] == sum(nodal_PTDF[jl, jn] * (g_scarcity[js,jn] - demand[js,jn] + ens_cm[js,jn]) for jn in JN))

        # mod.ext[:constraints][:cap_cm_global_balance] = @constraint(mod, [js in JS],
        #     0 == sum((g_scarcity[js,jn] - demand[js,jn] + ens_cm[js,jn]) for jn in JN))

    ############ Including exact projection constraints
    ## exchanges and flows
        for (jl,(i,j,_)) in enumerate(BRANCHES), js in JS
            if (zone_of[i],zone_of[j]) in TCONNECT
                mod.ext[:constraints][Symbol("ex_cm$(js)_$(jl)")] = @constraint(mod, flow_cm[js,jl] == ex_cm[js, (zone_of[i],zone_of[j])])
            elseif (zone_of[j],zone_of[i]) in TCONNECT
                mod.ext[:constraints][Symbol("ex_cm$(js)_$(jl)")] = @constraint(mod, flow_cm[js,jl] == -ex_cm[js, (zone_of[j],zone_of[i])])
            end
        end

        ## thermal limit
        mod.ext[:constraints][:thermal_limit] = @constraint(mod, [js in JS, jl in JL], -BRANCHES[jl][3] <= flow_cm[js,jl] <= BRANCHES[jl][3])

    elseif coupling == "ATC"     # each border is independently constrained by ATC 
        # ATC constraints

        # mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
        # sum(g_scarcity[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + sum(ens_cm[js,jn] for jn in JN if zone_of_idx[jn] == jz)
        # + cap_cm[jz] == sum(demand[js,jn] for jn in JN if zone_of_idx[jn] == jz))
           # # netposition = exports + imports (negative netposition --> export) 
        mod.ext[:constraints][:cap_cm_netposition] = @constraint(mod, [js in JS, jz in JZ],
            cap_cm[jz] ==
            sum(ex_cm[js,t] for t in TCONNECT if t[2] == zone_syms[jz])
            - sum(ex_cm[js,t] for t in TCONNECT if t[1] == zone_syms[jz])
            ) 
#  no redispatch in the ATC case
        mod.ext[:constraints][:cap_cm_atc_limit] = @constraint(mod, [js in JS, t in TCONNECT], 
            ATC[js][t][2] <= ex_cm[js,t] <= ATC[js][t][1])
    end

    return mod
end

