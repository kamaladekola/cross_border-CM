function solve_capacityIC_agent!(mod::Model)
    
    JZ = mod.ext[:sets][:JZ]
    JN = mod.ext[:sets][:JN]           # nodes
    JS = mod.ext[:sets][:JS]
    JL = mod.ext[:sets][:JL]           # lines

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar]
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]

    cap_cm = mod.ext[:variables][:cap_cm]
    ens_cm = mod.ext[:variables][:ens_cm]
    g_scarcity = mod.ext[:variables][:g_scarcity]



    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
    + sum(((ens_cm[js,jn])^2 * 1e10) for js in JS, jn in JN)
    )


    # for js in JS, jn in JN
    #     delete(mod, mod.ext[:constraints][:dispatch_feasibility][js,jn])
    # end
    # mod.ext[:constraints][:dispatch_feasibility] = @constraint(mod, [js in JS, jn in JN], 
    # g_scarcity[js, jn] ≤ CapCM_nodal[jn]) 


    optimize!(mod);

    return mod
end
