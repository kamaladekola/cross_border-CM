function solve_getATC!(mod::Model,
                       Cap_cm_nodal::AbstractVector,
                       d_scarcity::AbstractMatrix)

    m  = mod
    JN = mod.ext[:sets][:JN]
    JS = mod.ext[:sets][:JS]
    JV = mod.ext[:sets][:JV]
    JL = mod.ext[:sets][:JL]
    JZ = mod.ext[:sets][:JZ]
    zone_syms = mod.ext[:parameters][:zone_syms]
    nodes = mod.ext[:parameters][:nodes]
    nodal_PTDF = mod.ext[:parameters][:nodal_PTDF]
    CapCM_nodal = mod.ext[:parameters][:CapCM_nodal]
    zone_of = mod.ext[:parameters][:zone_of]
    d_scarcity  = mod.ext[:parameters][:d_scarcity] 



    TCONNECT = mod.ext[:parameters][:TCONNECT]
    atc_plus  = mod.ext[:variables][:atc_plus]
    atc_minus = mod.ext[:variables][:atc_minus]
    
    v_dayahead   = mod.ext[:variables][:v_dayahead]
    v_redispatch = mod.ext[:variables][:v_redispatch]
    f = mod.ext[:variables][:f]
    net_pos = mod.ext[:variables][:net_pos]
    epsilon = 1e-10

    zone_idx  = mod.ext[:parameters][:zone_of_idx]
    demand = mod.ext[:parameters][:getATC_demand]
    atc_results = Dict{Int,Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}}}()

    # mod.ext[:parameters][:CapCM_nodal] .= Cap_cm_nodal


    for js in JS
        demand .= d_scarcity[js, zone_idx] .* Cap_cm_nodal

        mod.ext[:objective] = @objective(m, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT)
        - epsilon * sum(atc_plus[t]^2 + atc_minus[t]^2 for t in TCONNECT)
        )


        for jv in JV, jl in JL
            delete(mod, mod.ext[:constraints][:getATC_nodal_balance][jv, jl])
        end
        mod.ext[:constraints][:getATC_nodal_balance] = @constraint(m, [jv in JV, jl in JL],
        # f[jv,jl] == sum(nodal_PTDF[jl, jn] * ((CapCM_nodal[jn] + ens[jn]) * (v_dayahead[jv,jn] + v_redispatch[jv,jn]) - demand[jn]) for jn in JN))
        f[jv,jl] == sum(nodal_PTDF[jl, jn] * ((CapCM_nodal[jn]) * (v_dayahead[jv,jn] + v_redispatch[jv,jn]) - demand[jn]) for jn in JN))

        for jv in JV, (jz, zsym) in enumerate(zone_syms)
            delete(mod, mod.ext[:constraints][:getATC_zonal_balance][jv, (jz, zsym)])
        end
        mod.ext[:constraints][:getATC_zonal_balance] = @constraint(m, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
            # sum((CapCM_nodal[jn] + ens[jn]) * v_dayahead[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym)
            sum((CapCM_nodal[jn]) * v_dayahead[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym)
            - net_pos[jv,jz] == sum(demand[jn] for jn in JN if zone_of[nodes[jn]] == zsym))


        for jv in JV, (jz, zsym) in enumerate(zone_syms)
            delete(mod, mod.ext[:constraints][:getATC_redispatch_limit][jv, (jz, zsym)])
        end
        mod.ext[:constraints][:getATC_redispatch_limit] =
            @constraint(m, [jv in JV, (jz, zsym) in enumerate(zone_syms)],
            # sum((CapCM_nodal[jn] + ens[jn]) * v_redispatch[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym) == 0)
            sum((CapCM_nodal[jn]) * v_redispatch[jv,jn] for jn in JN if zone_of[nodes[jn]] == zsym) == 0)

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