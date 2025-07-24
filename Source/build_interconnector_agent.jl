function build_interconnector_agent!(mod::Model)
    # Sets
    JH = mod.ext[:sets][:JH]           # hours
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    
    nodes = mod.ext[:parameters][:nodes] 
    W = mod.ext[:parameters][:w]             # hour weights
    # F_max = mod.ext[:parameters][:F_max]         # |L| thermal limits
    zonal_PTDF = mod.ext[:parameters][:zonal_PTDF]          # |L|×|Z|
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]      # |L|×|N| nodal PTDF
    λ_all = mod.ext[:parameters][:λ_all]         # |H|×|Z| EOM prices per zone
    # retrieve from ADMM_subroutine and sum all generation for each node.
    G_nodal = mod.ext[:parameters][:G_nodal]         # |H|×|N|
    D_nodal = mod.ext[:parameters][:D_nodal]         # |H|×|N|
    Y_nodal = mod.ext[:parameters][:Y_nodal]         # |H|×|N|

    g_bar_all = mod.ext[:parameters][:g_bar_all]     # |H|×|Z|
    ρ_all = mod.ext[:parameters][:ρ_all]         # |Z|
    zone_of = mod.ext[:parameters][:zone_of]       # Dict(node) → zone
    BRANCHES = mod.ext[:parameters][:BRANCHES]      #lines [(from,to,Fmax)]
    TCONNECT = mod.ext[:parameters][:TCONNECT]      # [(:A,:B),(:B,:C),…]
    zone_syms    = mod.ext[:parameters][:zone_syms] 
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]

    # net-position variable g[t,z] --> positive => import
    g = mod.ext[:variables][:g] = @variable(mod, g[jh=JH,jz=JZ], base_name = "netposition") # change to nodal aggregation
    g_red = mod.ext[:variables][:redispatch] = @variable(mod, g_red[jh in JH, jn in JN], base_name = "redispatch")
    flow  = mod.ext[:variables][:flow] = @variable(mod, flow[jh in JH, jl in JL], base_name = "flow")
    ex = mod.ext[:variables][:exchange] = @variable(mod, ex[jh in JH, t in TCONNECT], base_name = "exchange") # net exchange between zones

    ens_pos = mod.ext[:variables][:ens_pos] = @variable(mod, [jh in JH, jn in JN], lower_bound = 0)          # extra gen   (MW)
    ens_neg = mod.ext[:variables][:ens_neg] = @variable(mod, [jh in JH, jn in JN], lower_bound = 0)           # shed load   (MW)

    ens_IC = mod.ext[:expressions][:ens_IC] = @expression(mod, [jh in JH, jn in JN], ens_pos[jh,jn] - ens_neg[jh,jn])
    # ens_IC = mod.ext[:variables][:ens_IC] = @variable(mod, [jh in JH, jn in JN], base_name = "ens_IC")


    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    - sum(W[jh] * λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum(W[jh] * (ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 
    for jh in JH, jz in JZ))

    # Constraints --> exact projection onto the domain of FBMC to ensure network feasibility

    # crossborder exchanges: flow[jl] == e[(:A, :B)]  # import flows are stored as positive and export flows as negative
    for (jl,(i,j,_)) in enumerate(BRANCHES), jh in JH
        if (zone_of[i],zone_of[j]) in TCONNECT
            mod.ext[:constraints][Symbol("ex_$(jh)_$(jl)")] =  @constraint(mod, flow[jh,jl] ==  ex[jh, (zone_of[i],zone_of[j])])
        elseif (zone_of[j],zone_of[i]) in TCONNECT
            mod.ext[:constraints][Symbol("ex_$(jh)_$(jl)")] = @constraint(mod, flow[jh,jl] == -ex[jh, (zone_of[j],zone_of[i])])
        end
    end

    # thermal limit constraint
    mod.ext[:constraints][:thermal_limit] = @constraint(mod, [jh in JH, jl in JL], -BRANCHES[jl][3] ≤ flow[jh,jl] ≤ BRANCHES[jl][3])

    # netposition = exports + imports (negative netposition --> export) 
    mod.ext[:constraints][:net_pos] = @constraint(mod, [jh in JH, jz in JZ],
      g[jh,jz] ==
        sum(ex[jh, t] for t in TCONNECT if t[2] == zone_syms[jz])
        - sum(ex[jh, t] for t in TCONNECT if t[1] == zone_syms[jz])
    )

    # mod.ext[:constraints][:net_pos] = @constraint(mod, [jh in JH, jz in JZ], g[jh,jz] == 0) # killing interconnector


    # nodal balance constraint: flow = generation + redispatch - demand --> for all nodes
    mod.ext[:constraints][:nodal_balance] = @constraint(mod, [jh in JH, jl in JL],
        flow[jh,jl] ==
          sum(nodal_PTDF[jl, jn] * (G_nodal[jh,jn] + g_red[jh,jn] - D_nodal[jh,jn] + ens_IC[jh,jn]) for jn in JN)
    )

    mod.ext[:constraints][:np_balance] = @constraint(mod, [jh in JH], sum(g[jh,jz] for jz in JZ) == 0) # redundant (saveguard)

    # generation + redispatch <= installed capacity (Y_nodal) * availability factor -> this matters especially for renewables
    # redispatch constraints
    mod.ext[:constraints][:redispatch_limit] = @constraint(mod, [jh in JH, jn in JN], 
    0 <= G_nodal[jh,jn] + g_red[jh,jn] + ens_IC[jh,jn] <= Y_nodal[jh,jn]) 

    mod.ext[:constraints][:redispatch_balance] = @constraint(mod, [jh in JH, jz in JZ],
      sum(g_red[jh, jn] + ens_IC[jh,jn] for jn in JN if zone_of_idx[jn] == jz) == 0)

    # mod.ext[:constraints][:redispatch_balance] = @constraint(mod, [jh in JH, jz in JZ],
    #   -1 <= sum(g_red[jh, jn] + ens_IC[jh,jn] for jn in JN if zone_of_idx[jn] == jz) <= 1)
    # bound on redispatch not needed since it is already bounded by the capacity limit constraint


    return mod
end
