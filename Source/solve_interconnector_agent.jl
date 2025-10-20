function solve_interconnector_agent!(mod::Model)

    # Sets
    JH = mod.ext[:sets][:JH]           # hours
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes

    # Parameters
    W = mod.ext[:parameters][:w]                            # hour weights
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]         # |L|×|N| nodal PTDF
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
    g = mod.ext[:variables][:g]                    # net position: positive => import
    r = mod.ext[:variables][:r]                    # nodal injections
    flow  = mod.ext[:variables][:flow]             # line flows

    g_bar = mod.ext[:variables][:g_bar]             # nodal generation
    y_bar = mod.ext[:variables][:y_bar]             # nodal capacity allocation
    s = mod.ext[:variables][:s]

    # Objective
    mod.ext[:objective] = @objective(mod, Min, 
    - sum(W[jh] * λ_all[jh,jz] * g[jh,jz] for jh in JH, jz in JZ)
    + sum(W[jh] * (ρ_all[jz]/2) * (g[jh,jz] - g_bar_all[jh,jz])^2 for jh in JH, jz in JZ)
    + sum(rc * s[jn] for jn in JN)
    )

    # -------------------------------------------------------------
    # redefining constraints with updated parameters
    # --------------------------------------------------------------

    # nodal balance constraint: flow = generation - demand
    for jh in JH, jn in JN
        delete(mod, mod.ext[:constraints][:nodal_balance][jh,jn])
    end
    mod.ext[:constraints][:nodal_balance] = @constraint(mod, [jh in JH, jn in JN], 
        r[jh, jn] == g_bar[jh, jn] - D_nodal[jh,jn])

    # nodal capacity limit
    for jh in JH, jn in JN
        delete(mod, mod.ext[:constraints][:cap_limit][jh,jn])
    end
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh in JH, jn in JN], (g_bar[jh, jn] - Y_nodal[jh, jn]) <= y_bar[jn] + s[jn]) # to be checked

    # zonal capacity allocation to nodes
    for jz in JZ
        delete(mod, mod.ext[:constraints][:capacity_allocation][jz])
    end
    mod.ext[:constraints][:capacity_allocation] = @constraint(mod, [jz in JZ], Y_zonal[jz] == sum(y_bar[jn] for jn in JN if zone_of_idx[jn] == jz))

    optimize!(mod);

    return mod
end
