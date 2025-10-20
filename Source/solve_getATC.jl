function solve_getATC!(mod::Model)

    JN = mod.ext[:sets][:JN]
    JS = mod.ext[:sets][:JS]
    JV = mod.ext[:sets][:JV]
    JL = mod.ext[:sets][:JL]
    JZ = mod.ext[:sets][:JZ]
    zone_syms = mod.ext[:parameters][:zone_syms]
    nodes = mod.ext[:parameters][:nodes]
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    zone_of_node = mod.ext[:parameters][:zone_of_node]
    d_scarcity  = mod.ext[:parameters][:d_scarcity] 
    Cap_Demand_nodal = mod.ext[:parameters][:Cap_Demand_nodal]
    reserve_cost = mod.ext[:parameters][:reserve_cost]

    epsilon = 1e-8

    TCONNECT = mod.ext[:parameters][:TCONNECT]
    atc_plus  = mod.ext[:variables][:atc_plus]
    atc_minus = mod.ext[:variables][:atc_minus]
    
    v_dayahead   = mod.ext[:variables][:v_dayahead]
    v_redispatch = mod.ext[:variables][:v_redispatch]
    f = mod.ext[:variables][:f]
    net_pos = mod.ext[:variables][:net_pos]
    reserve = mod.ext[:variables][:network_reserve]

    atc_results = Dict{Int,Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}}}()


    for js in JS

        demand = mod.ext[:expressions][:v_demand] = @expression(mod, [jn in JN],
            d_scarcity[js, jn] * Cap_Demand_nodal[jn])

        mod.ext[:objective] = @objective(mod, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT)
        - epsilon * sum(atc_plus[t]^2 + atc_minus[t]^2 for t in TCONNECT)
        - reserve_cost * sum(reserve[jn] for jn in JN)
        )


        for jv in JV, jl in JL
            delete(mod, mod.ext[:constraints][:getATC_nodal_balance][jv, jl])
        end
        mod.ext[:constraints][:getATC_nodal_balance] = @constraint(mod, [jv in JV, jl in JL],
        f[jv,jl] == sum(nodal_PTDF[jl, jn] * (CapCM_nodal[jn] * (v_dayahead[jv,jn] + v_redispatch[jv,jn]) + reserve[jn] - demand[jn]) for jn in JN))

        for jv in JV, (jz, zsym) in enumerate(zone_syms)
            delete(mod, mod.ext[:constraints][:getATC_zonal_balance][jv, (jz, zsym)])
        end

        mod.ext[:constraints][:getATC_zonal_balance] = @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
        sum(CapCM_nodal[jn] * v_dayahead[jv,jn] for jn in JN if zone_of_node[nodes[jn]] == zsym)
        - net_pos[jv,jz] == sum(demand[jn] - reserve[jn]  for jn in JN if zone_of_node[nodes[jn]] == zsym)
        )

        for jv in JV, (jz, zsym) in enumerate(zone_syms)
            delete(mod, mod.ext[:constraints][:getATC_redispatch_limit][jv, (jz, zsym)])
        end
        mod.ext[:constraints][:getATC_redispatch_limit] =
        @constraint(mod, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
        sum((CapCM_nodal[jn]) * v_redispatch[jv,jn] for jn in JN if zone_of_node[nodes[jn]] == zsym) == 0)

        optimize!(mod)

        status = JuMP.termination_status(mod)
        if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
            error("ATC model failed in scenario $js – status = $status")
        end



        plus   = value.(atc_plus)
        minus  = value.(atc_minus)

        atc_results[js] = Dict(t => (abs(plus[t]), -abs(minus[t])) for t in TCONNECT)
    end
    return atc_results
end