function build_interconnector_agent!(mod::Model)

    # Sets
    JH = mod.ext[:sets][:JH]           # hours
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes

    # Parameters
    W = mod.ext[:parameters][:w]                            # hour weights
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]          # |L|×|N| nodal PTDF
    λ_all = mod.ext[:parameters][:λ_all]                    # |H|×|Z| EOM prices per zone

    # retrieved from ADMM_subroutine
    D_nodal = mod.ext[:parameters][:D_nodal]                # |H|×|N|
    Y_nodal = mod.ext[:parameters][:Y_nodal]                # |H|×|N|
    Y_zonal = mod.ext[:parameters][:Y_zonal]                # |Z| zonal capacity allocation
    g_bar_all = mod.ext[:parameters][:g_bar_all]            # |H|×|Z|
    ρ_all = mod.ext[:parameters][:ρ_all]                    # |Z|
    BRANCHES = mod.ext[:parameters][:BRANCHES]              #lines [(from,to,Fmax)]
    Fmax = [BRANCHES[jl][3] for jl in JL]
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]
    rc = mod.ext[:parameters][:reserve_cost]

    # Variables
    g = mod.ext[:variables][:g] = @variable(mod, g[jh=JH,jz=JZ], base_name = "netposition")                     # net position: positive => import
    r = mod.ext[:variables][:r] = @variable(mod, r[jh in JH, jn in JN], base_name = "nodal_injection")          # nodal injections
    flow  = mod.ext[:variables][:flow] = @variable(mod, flow[jh in JH, jl in JL], base_name = "flow")           # line flows
    
    g_bar = mod.ext[:variables][:g_bar] = @variable(mod, [jh=JH, jn=JN], lower_bound=0, base_name="nodal_generation")          # nodal generation
    y_bar = mod.ext[:variables][:y_bar] = @variable(mod, [jn=JN], lower_bound=0, base_name="nodal_capacity")    # nodal capacity allocation
    s = mod.ext[:variables][:s] = @variable(mod, [jn in JN], lower_bound = 0, base_name = "network_reserve")    # reserve provision to ensure feasibility

    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    - sum(W[jh] * λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum(W[jh] * (ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 for jh in JH, jz in JZ)
    + sum(rc * s[jn] for jn in JN)
    )

    # -------------------------------------------------------------
    # Constraints --> FBMC using exact projection
    # --------------------------------------------------------------

    # nodal balance constraint: flow = generation - demand
    mod.ext[:constraints][:nodal_balance] = @constraint(mod, [jh in JH, jn in JN], 
        r[jh, jn] == g_bar[jh, jn] - D_nodal[jh,jn])

    # nodal capacity limit
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh in JH, jn in JN], (g_bar[jh, jn] - Y_nodal[jh, jn]) <= y_bar[jn] + s[jn]) # to be checked

    # zonal capacity allocation to nodes
    mod.ext[:constraints][:capacity_allocation] = @constraint(mod, [jz in JZ], Y_zonal[jz] == sum(y_bar[jn] for jn in JN if zone_of_idx[jn] == jz))

    # Aggregate nodal flows to zonal net positions
    mod.ext[:constraints][:net_pos] = @constraint(mod, [jh in JH, jz in JZ],
        g[jh, jz] == - sum(r[jh, jn] for jn in JN if zone_of_idx[jn] == jz)
    )

    ### killing interconnector
    # mod.ext[:constraints][:net_pos] = @constraint(mod, [jh in JH, jz in JZ], g[jh,jz] == 0)

    # System balance (sum injections = 0)
    mod.ext[:constraints][:sys_balance] = @constraint(mod, [jh in JH],
        sum(r[jh, jn] for jn in JN) == 0
    )

    # DC power flow mapping: flows = PTDF * injections
    mod.ext[:constraints][:flows] = @constraint(mod, [jh in JH, jl in JL],
        flow[jh, jl] == sum(nodal_PTDF[jl, jn] * r[jh, jn] for jn in JN)
    )

    # thermal limits
    mod.ext[:constraints][:thermal] = @constraint(mod, [jh in JH, jl in JL],
        -Fmax[jl] <= flow[jh, jl] <= Fmax[jl]
    )

    return mod
end
