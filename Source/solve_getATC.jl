function solve_getATC!(mod::Model; hour::Union{Int,Nothing}=nothing)

        # Sets
    JH = mod.ext[:sets][:JH]           # hours
    JZ = mod.ext[:sets][:JZ]           # zones
    JL = mod.ext[:sets][:JL]           # lines
    JN = mod.ext[:sets][:JN]           # nodes
    JV = mod.ext[:sets][:JV]             # vertices
    
    nodes = mod.ext[:parameters][:nodes] 
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    D_nodal = mod.ext[:parameters][:D_nodal] 
    Y_nodal = mod.ext[:parameters][:Y_nodal]


    zone_of = mod.ext[:parameters][:zone_of]
    BRANCHES = mod.ext[:parameters][:BRANCHES]
    TCONNECT = mod.ext[:parameters][:TCONNECT]
    zone_syms    = mod.ext[:parameters][:zone_syms] 
    zone_of_idx = mod.ext[:parameters][:zone_of_idx]

    atc_plus = mod.ext[:variables][:atc_plus]
    atc_minus = mod.ext[:variables][:atc_minus]


    v_dayahead = mod.ext[:variables][:v_dayahead]
    v_redispatch = mod.ext[:variables][:v_redispatch]
    net_pos = mod.ext[:variables][:net_pos]
    e = mod.ext[:variables][:e]
    
    @objective(m, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT))

    for js in JS, jn in JN
        delete(mod, mod.ext[:constraints][:dispatch_feasibility][js,jn])
    end
    mod.ext[:constraints][:dispatch_feasibility] = @constraint(mod, [js in JS, jn in JN], 
    g_scarcity[js, jn] ≤ CapCM_nodal[jn]) 

        for jz in JZ
            for jv in JV, jn in JN
                delete(mod, mod.ext[:constraints][:getATC_zonal_balance][jv,jn])
                delete(mod, mod.ext[:constraints][:getATC_redispatch_limit][jv,jn])
            end

            mod.ext[:constraints][:getATC_zonal_balance] = @constraint(m, sum(Y_nodal[jn] *v_dayahead[jv,jn] for jn in JN if zone_idx_of[jn] == jz) - net_pos[v,jz]
                == sum(D_nodal[jn] for jn in nodes if zone_idx_of[jn]==jz))
                           
            mod.ext[:constraints][:getATC_redispatch_limit] = @constraint(m, sum(Y_nodal[jn] * v_redispatch[v,jn] for jn in JN if zone_idx_of[jn] == jz) == 0)                        
        end

    optimize!(mod)


    atc_plus  = value.(mod.ext[:variables][:atc_plus])
    atc_minus = value.(mod.ext[:variables][:atc_minus])

    atc_vector = [(abs(atc_plus[t]), -abs(atc_minus[t])) for t in mod.ext[:parameters][:TCONNECT]]

    return mod, atc_vector
end