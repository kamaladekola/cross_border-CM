function build_getATC!(mod::Model)

    JV = mod.ext[:sets][:JV]
    JH = mod.ext[:sets][:JH]
    JL = mod.ext[:sets][:JL]
    JN = mod.ext[:sets][:JN]
    JS = mod.ext[:sets][:JS]

    TCONNECT    = mod.ext[:parameters][:TCONNECT]
    signs       = mod.ext[:parameters][:signs]
    nodes       = mod.ext[:parameters][:nodes]
    BRANCHES    = mod.ext[:parameters][:BRANCHES]
    zone_syms   = mod.ext[:parameters][:zone_syms]
    zone_of     = mod.ext[:parameters][:zone_of]
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]

    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    nodal_PTDF  = mod.ext[:parameters][:nodal_PTDF]
    d_scarcity  = mod.ext[:parameters][:d_scarcity] 

    demand = mod.ext[:parameters][:Cap_Demand_nodal]

    epsilon = 1e-8
    M = 1e10

    
    s_pos = mod.ext[:variables][:s_pos] = @variable(mod, [jn in JN], lower_bound = 0)          # extra gen   (MW)
    s_neg = mod.ext[:variables][:s_neg] = @variable(mod, [jn in JN], lower_bound = 0)           # shed load   (MW)

    ens = mod.ext[:expressions][:ens_atc] = @expression(mod, [jn=JN], s_pos[jn] - s_neg[jn])


    atc_plus  = mod.ext[:variables][:atc_plus]  =
        @variable(mod, [t in TCONNECT], base_name = "atc_plus")

    atc_minus = mod.ext[:variables][:atc_minus] =
        @variable(mod, [t in TCONNECT], base_name = "atc_minus")

    v_dayahead   = mod.ext[:variables][:v_dayahead]   =
        @variable(mod, [jv in JV, jn in JN], base_name = "v_dayahead", lower_bound = 0, upper_bound = 1)

    v_redispatch = mod.ext[:variables][:v_redispatch] =
        @variable(mod, [jv in JV, jn in JN], base_name = "v_redispatch", lower_bound = -1, upper_bound = 1)

    f = mod.ext[:variables][:f] = @variable(mod, [jv in JV, jl in JL], base_name = "f")

    net_pos = mod.ext[:variables][:net_pos] = @variable(mod, [jv in JV, jz in 1:length(zone_syms)], base_name = "net_pos")

    e = mod.ext[:variables][:e] = @variable(mod, [jv in JV, t in TCONNECT], base_name = "e")



    # mod.ext[:objective] = @objective(mod, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT) - sum(ens[jn] for jn in JN)^2 * 1e10)
    mod.ext[:objective] = @objective(mod, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT)
        - epsilon * sum(atc_plus[t]^2 + atc_minus[t]^2 for t in TCONNECT)
        - M * sum((s_pos[jn] + s_neg[jn]) for jn in JN)
        )

    # @constraint(mod, [t in TCONNECT], -atc_minus[t] <= atc_plus[t])

    # CONSTRAINTS
    @constraint(mod, [jv in JV, jn in JN], v_dayahead[jv, jn] + v_redispatch[jv, jn] >= 0)
    @constraint(mod, [jv in JV, jn in JN], v_dayahead[jv, jn] + v_redispatch[jv, jn] <= 1)

    # map vertex signs to the +/- ATC variables
    @constraint(mod, [jv in JV, (k, t) in enumerate(TCONNECT)],
        e[jv, t] == (signs[jv][k] == 1 ?  atc_plus[t] : -atc_minus[t]))

    # thermal limits
    @constraint(mod, [jv in JV, jl in JL],
        -BRANCHES[jl][3] <= f[jv,jl] <= BRANCHES[jl][3])

    mod.ext[:constraints][:getATC_nodal_balance] = @constraint(mod, [jv in JV, jl in JL],
        # f[jv,jl] == sum(nodal_PTDF[jl, jn] * ((CapCM_nodal[jn]) * (v_dayahead[jv,jn] + v_redispatch[jv,jn]) - demand[jn]) for jn in JN))
        f[jv,jl] == sum(nodal_PTDF[jl, jn] * (CapCM_nodal[jn] * (v_dayahead[jv,jn] + v_redispatch[jv,jn]) - demand[jn]  + ens[jn]) for jn in JN))

    @constraint(mod, [jv in JV, t in TCONNECT],
        e[jv,t] == sum(f[jv,jl] * ((zone_of[BRANCHES[jl][1]], zone_of[BRANCHES[jl][2]]) == t  ?  1 :
                 (zone_of[BRANCHES[jl][2]], zone_of[BRANCHES[jl][1]]) == t  ? -1 : 0) for jl in JL))

    # zone net positions
    @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
        net_pos[jv,jz] ==
            sum(e[jv,t] for t in TCONNECT if t[1] == zsym) -
            sum(e[jv,t] for t in TCONNECT if t[2] == zsym))

    mod.ext[:constraints][:getATC_zonal_balance] = @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
        sum(CapCM_nodal[jn] * v_dayahead[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym)
        # - net_pos[jv,jz] == sum(demand[jn] for jn in JN if zone_of[nodes[jn]] == zsym)
        - net_pos[jv,jz] == sum(demand[jn] - ens[jn] for jn in JN if zone_of[nodes[jn]] == zsym))

    mod.ext[:constraints][:getATC_redispatch_limit] =
        @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
        # sum((CapCM_nodal[jn] + ens[jn]) * v_redispatch[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym) == 0)
        sum((CapCM_nodal[jn]) * v_redispatch[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym) == 0)


    return mod
end