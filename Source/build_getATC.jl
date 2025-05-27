
function build_getATC!(mod::Model, hour=20)

    # Sets
    JH = mod.ext[:sets][:JH]           # hours
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    JV = mod.ext[:sets][:JV]             # vertices
    
    TCONNECT   = mod.ext[:parameters][:TCONNECT]
    signs      = mod.ext[:parameters][:signs]
    nodes      = mod.ext[:parameters][:nodes]
    D_nodal    = mod.ext[:parameters][:D_nodal] # use demand at scarcity scenarios - d_scarcity
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    zone_idx_of = mod.ext[:parameters][:zone_idx_of]
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    zone_of    = mod.ext[:parameters][:zone_of]
    BRANCHES   = mod.ext[:parameters][:BRANCHES]
    zone_syms  = mod.ext[:parameters][:zone_syms]

    # ATC variables - two for  each interconnector
    # solve and obtain atc for each scenario
    atc_plus = mod.ext[:variables][:atc_plus] = @variable(m, atc_plus[t in TCONNECT], base_name = "atc_plus")
    atc_minus = mod.ext[:variables][:atc_minus] = @variable(m, atc_minus[t in TCONNECT], base_name = "atc_minus")


    v_dayahead = mod.ext[:variables][:v_dayahead] = @variable(m, 0 <= v_dayahead[jv in JV, jn in JN] <= 1, base_name = "v_dayahead")
    v_redispatch = mod.ext[:variables][:v_redispatch] = @variable(m, 0 <= v_redispatch[jv in JV, jn in JN] <= 1, base_name = "v_redispatch")
    f = mod.ext[:variables][:f] = @variable(m, f[jv in JV, jl in JL] >= 0, base_name = "f")
    net_pos = mod.ext[:variables][:net_pos] = @variable(m, net_pos[jv in JV, jz in JZ] >= 0, base_name = "net_pos")
    e = mod.ext[:variables][:e] = @variable(m, e[jv in JV, t in TCONNECT] >= 0, base_name = "e")
    
    ###########################################################################################################
    demand = mod.ext[:expressions][:demand] = @expression(mod, [js in JS, jn in JN],
            d_scarcity[js, zone_of_idx[jn]] * CapCM_nodal[jn]) 
            
    # to do: define nodal demand as a fraction of total system capacity
    ###########################################################################################################

    # Objective
    @objective(m, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT))
    
    # Constraints
    @constraint(mod, [jv in JV, jn in JN],
      v_dayahead[jv,jn] + v_redispatch[jv,jn] >= 0
    )
    @constraint(mod, [jv in JV, jn in JN],
      v_dayahead[jv,jn] + v_redispatch[jv,jn] <= 1
    )


    # For each vertex
    # for (v,s) in enumerate(signs)
        # @constraint(m, [t in TCONNECT], e[v,t] == -atc_minus[t])
        # @constraint(m, [t in TCONNECT], e[v,t] == atc_plus[t])

    @constraint(mod, [jv in JV, (k,t) in enumerate(TCONNECT)],
      e[jv,t] == ( signs[jv][k] == 1 ? atc_plus[t] : -atc_minus[t] )
    )

    @constraint(mod, [jv in JV, jl in JL],
          -BRANCHES[jl][3] <= f[jv,jl] <= BRANCHES[jl][3]
    )

    mod.ext[:constraints][:getATC_nodal_balance] = @constraint(mod, [jv in JV, jl in JL],
      f[jv,jl] == sum(nodal_PTDF[jl, jn] * ( Y_nodal[jn] * (v_dayahead[jv,jn] + v_redispatch[jv,jn])
              - D_nodal[jn]) for jn in JN))


    @constraint(mod, [jv in JV, t in TCONNECT],
      e[jv,t] == sum(f[jv,jl] * ((zone_of[BRANCHES[jl][1]], zone_of[BRANCHES[jl][2]]) == t ? 1
                    : (zone_of[BRANCHES[jl][2]], zone_of[BRANCHES[jl][1]]) == t ? -1
                    : 0) for jl in JL))

    @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
      net_pos[jv,jz] == 
      sum(e[jv,t] for t in TCONNECT if t[1] == zsym)
      - sum(e[jv,t] for t in TCONNECT if t[2] == zsym)
      )


    mod.ext[:constraints][:getATC_zonal_balance] = @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
      sum(Y_nodal[jn] * v_dayahead[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym)
      - net_pos[jv,jz] == sum(D_nodal[jn] for jn in JN if zone_of[nodes[jn]] == zsym)
    )

    mod.ext[:constraints][:getATC_redispatch_limit] =  @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
      sum(Y_nodal[jn] * v_redispatch[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym) == 0
    )


    return mod
end
# end

 

