function build_capacityIC_agent!(mod::Model)
    # Sets
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]           # scarcity scenarios
    
    # Parameters
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    zone_syms = mod.ext[:parameters][:zone_syms]
    node_syms = mod.ext[:parameters][:nodes]
    zone_idx = Dict(z => i for (i, z) in enumerate(zone_syms))
    BRANCHES = mod.ext[:parameters][:BRANCHES]
    Fmax = [BRANCHES[jl][3] for jl in JL]
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    scarcity_matrix = mod.ext[:parameters][:scarcity_matrix]
    CapCM_zonal = mod.ext[:parameters][:CapCM_zonal]
    Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal]
    rc = mod.ext[:parameters][:reserve_cost]
    coupling = mod.ext[:parameters][:coupling]

    zone_of_node = mod.ext[:parameters][:zone_of_node]
    
    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar]
    λ_CM = mod.ext[:parameters][:λ_CM]
    ρ_CM = mod.ext[:parameters][:ρ_CM]
    
    y_bar_nodal = mod.ext[:parameters][:y_bar_nodal]
    called = mod.ext[:parameters][:called]
    
    # Variables
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], base_name = "netposition")                     # net position: positive => import
    r_cm = mod.ext[:variables][:r_cm] = @variable(mod, [js=JS, jn=JN], base_name = "nodal_injection")          # nodal injections
    flow_cm  = mod.ext[:variables][:flow_cm] = @variable(mod, [js=JS, jl=JL], base_name = "flow")           # line flows
    s_cm = mod.ext[:variables][:s_cm] = @variable(mod, [js=JS, jn=JN], lower_bound = 0, base_name = "network_reserve")    # reserve provision to ensure feasibility
    CapCM_nodal = mod.ext[:variables][:CapCM_nodal] = @variable(mod, [jn=JN], lower_bound=0, base_name="nodal_capacity")    # nodal capacity allocation
    g_scar = mod.ext[:variables][:g_scar] = @variable(mod, [js=JS, jn=JN], lower_bound=0, base_name="nodal_generation")          # nodal generation


    # Expressions
    # demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
    #     d_scarcity[js, jn] * Cap_Demand_nodal[jn])
    
    demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
        scarcity_matrix[js, zone_idx[zone_of_node[node_syms[jn]]]] * Cap_Demand_nodal[jn])
        
    # Objective function
    mod.ext[:objective] = @objective(mod, Min,
        - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
        + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
        + sum(rc * s_cm[js, jn] for js in JS, jn in JN)
        )

    
    # -------------------
    # Constraints
    # -------------------
    # g_scar ≤ CapCM_nodal + s_res
    mod.ext[:constraints][:cap_limit] =
        @constraint(mod, [js in JS, jn in JN], g_scar[js, jn] <= CapCM_nodal[jn])

    # CapCM_nodal <= y_bar_nodal
    mod.ext[:constraints][:capcm_limit] = @constraint(mod, [jn in JN], CapCM_nodal[jn] <= y_bar_nodal[jn])

    # zonal capacity allocation to nodes
    mod.ext[:constraints][:capacity_allocation] = @constraint(mod, [jz in JZ],
        sum(CapCM_nodal[jn] for jn in JN if zone_of_idx[jn] == jz) == CapCM_zonal[jz])

    if coupling == "FB"

        # nodal balance constraint: for all scenarios, nodal injection = generation - demand
        mod.ext[:constraints][:nodal_balance] = @constraint(mod, [js in JS, jn in JN], 
            r_cm[js, jn] == g_scar[js, jn] - (demand[js,jn] - s_cm[js, jn]))

        # Zonal net positions from nodal injections (import-positive)
        #    cap_cm[z] >= - Σ_{n∈z} r_cm[js,n]  for ALL scenarios js  (robust deliverability)

        # mod.ext[:constraints][:netpos] =
        #     @constraint(mod, [js in JS, jz in JZ],
        #         cap_cm[jz] >= - sum(r_cm[js, jn] for jn in JN if zone_of_idx[jn] == jz)
        #     )

        np_s = mod.ext[:expressions][:np_s] = @expression(mod, [js in JS, jz in JZ], -sum(r_cm[js, jn] for jn in JN if zone_of_idx[jn] == jz))

        sys_called = Dict(js => any(called[js, jz] for jz in JZ) for js in JS)

        mod.ext[:constraints][:netpos_imp] = @constraint(mod, [js in JS, jz in JZ; called[js, jz]],
            np_s[js, jz] >= cap_cm[jz])

        mod.ext[:constraints][:netpos_exp] = @constraint(mod, [js in JS, jz in JZ; sys_called[js] && !called[js, jz]],
            np_s[js, jz] <= cap_cm[jz])

        # deliverability floor + demand shortfall for called zones
        mod.ext[:constraints][:netpos] = @constraint(mod, [js in JS, jz in JZ; called[js, jz]],
            np_s[js, jz] >= cap_cm[jz] + sum(demand[js, jn] - Cap_Demand_nodal[jn] for jn in JN if zone_of_idx[jn] == jz))

        # System balance (sum injections = 0) for all scenarios
        mod.ext[:constraints][:sys_bal] =
            @constraint(mod, [js in JS], sum(r_cm[js, jn] for jn in JN) == 0)

        # Overall net position sums to zero
        mod.ext[:constraints][:global_bal] =
            @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0)

        # flows = PTDF * nodal injections
        mod.ext[:constraints][:flows] =
            @constraint(mod, [js in JS, jl in JL], flow_cm[js, jl] == sum(nodal_PTDF[jl, jn] * r_cm[js, jn] for jn in JN))

        # Thermal limits
        mod.ext[:constraints][:thermal] =
            @constraint(mod, [js in JS, jl in JL], -Fmax[jl] <= flow_cm[js, jl] <= Fmax[jl])

    elseif coupling == "ATC"
        TCONNECT = mod.ext[:parameters][:TCONNECT]
        ATC = mod.ext[:parameters][:ATC]

        ex_cm = mod.ext[:variables][:ex_cm] = @variable(mod, [js in JS, t in TCONNECT], base_name = "ex_cm")

        # scenario net position, import-positive; Σ_z np_s[js,:] == 0 identically (= sys_bal)
        np_s = mod.ext[:expressions][:np_s] = @expression(mod, [js in JS, jz in JZ],
            sum(ex_cm[js, t] for t in TCONNECT if t[2] == zone_syms[jz]) -
            sum(ex_cm[js, t] for t in TCONNECT if t[1] == zone_syms[jz]))

        # zonal balance = FB nodal_balance summed over the zone
        mod.ext[:constraints][:zonal_balance] = @constraint(mod, [js in JS, jz in JZ],
            np_s[js, jz] == sum(demand[js, jn] - s_cm[js, jn] - g_scar[js, jn] for jn in JN if zone_of_idx[jn] == jz))

        mod.ext[:constraints][:global_bal] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0)

        sys_called = Dict(js => any(called[js, jz] for jz in JZ) for js in JS)

        # called zones: deliverability floor
        mod.ext[:constraints][:netpos_imp] = @constraint(mod, [js in JS, jz in JZ; called[js, jz]],
            np_s[js, jz] >= cap_cm[jz])

        # uncalled zones in a scenario with a call: honour the obligation
        mod.ext[:constraints][:netpos_exp] = @constraint(mod, [js in JS, jz in JZ; sys_called[js] && !called[js, jz]],
            np_s[js, jz] <= cap_cm[jz])

        # called zones: deliverability floor + demand shortfall
        mod.ext[:constraints][:netpos_imp] = @constraint(mod, [js in JS, jz in JZ; called[js, jz]],
            np_s[js, jz] >= cap_cm[jz] + sum(demand[js, jn] - Cap_Demand_nodal[jn] for jn in JN if zone_of_idx[jn] == jz))

        # ATC limits per scenario
        mod.ext[:constraints][:cap_cm_atc_limit] = @constraint(mod, [js in JS, t in TCONNECT],
            ATC[t][2] <= ex_cm[js, t] <= ATC[t][1])
        else
        # constrain cap_cm to zero
        @constraint(mod, [jz in JZ], cap_cm[jz] == 0)
    end
    
    return mod
end