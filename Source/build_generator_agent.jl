function build_generator_agent!(mod::Model, m::String, zones::Vector{String})

    zone, _ = parse_agent_name(m)
    home_zone = findfirst(isequal(zone), zones)
    # Extract sets
    JH = mod.ext[:sets][:JH]   # timesteps
    JZ = mod.ext[:sets][:JZ]   # zones
    JN = mod.ext[:sets][:JN]   # nodes

   # Extract parameters
    AF = mod.ext[:timeseries][:AF]
    W = mod.ext[:parameters][:w]
    zone_of_idx   = mod.ext[:parameters][:zone_of_idx]
    zone_syms = mod.ext[:parameters][:zone_syms]
    node_share = mod.ext[:parameters][:node_share]
    nodes= mod.ext[:parameters][:nodes]



    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM]                # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]                # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]                # rho-value in ADMM related to EOM auctions
    cap_bar = mod.ext[:parameters][:cap_bar]            # element in ADMM penalty term related to capacity markets
    ρ_CM = mod.ext[:parameters][:ρ_CM]
    λ_CM = mod.ext[:parameters][:λ_CM]
    σ_CM = mod.ext[:parameters][:σ_CM]                  # 1 if capacity markets are active, 0 otherwise
    PM = mod.ext[:parameters][:participation_matrix]
    DF = mod.ext[:parameters][:derating_factor]

    # Investment parameters
    I = mod.ext[:parameters][:I]
    max_cap  = mod.ext[:parameters][:max_cap]
    y_init   = mod.ext[:parameters][:C]   # existing capacity

    cap_cm_nodal = mod.ext[:variables][:cap_cm_nodal] = @variable(mod, [jn in JN], lower_bound=0, base_name = "cap_cm_nodal")
    g_nodal = mod.ext[:variables][:g_nodal] = @variable(mod, [jh in JH, jn in JN], lower_bound=0, base_name = "g_nodal")
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], lower_bound=0, base_name="generation")                                # generation
    y = mod.ext[:variables][:y] = @variable(mod, lower_bound=y_init, upper_bound=max_cap, base_name="InstalledCapacity")       # installed capacity
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], lower_bound = 0, base_name = "capacity_offered")                     # capacity offered in capacity markets
    


    # Objective => minimize GenCo costs
    mod.ext[:objective] = @objective(mod, Min,
        # + sum(W[jh] * A/2*g[jh]^2 for jh in JH)                       # cost function for generation with weights
        + sum(W[jh] * B*g[jh] for jh in JH)
        - sum(W[jh] * λ_EOM[jh]*g[jh] for jh in JH)                   # revenue from EOM with weights
        + I * (y - y_init)                                            # investment cost
        # - σ_CM * sum(λ_CM[jz] * cap_cm[jz] * PM[m][zones[jz]] * DF[m][zones[jz]] for jz in JZ)        
        - σ_CM * sum(λ_CM[home_zone] * cap_cm[jz] * PM[m][zones[jz]] * DF[m][zones[jz]] for jz in JZ) # when all capacity in a zone is remunerated at the same price
        + sum(W[jh] * ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)     # ADMM penalty term for EOM clearing with weights
        # + σ_CM * sum(ρ_CM[jz]/2 * PM[m][zones[jz]] * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ) # ADMM penalty term for capacity markets
        + σ_CM * ρ_CM[home_zone]/2 * (sum(cap_cm[jz] * PM[m][zones[jz]] for jz in JZ) - cap_bar[home_zone])^2 # ADMM penalty term for CM when all capacity in a zone is remunerated at the same price
    )

    # # aggregate nodal to zone level
    mod.ext[:constraints][:g_link] = @constraint(mod, [jh in JH, jn in JN],
    g_nodal[jh,jn] == get(node_share, nodes[jn], 0.0) * g[jh])

    # # available capacity in node
    y_nodal = mod.ext[:expressions][:y_nodal] = @expression(mod, [jh in JH, jn in JN],
            AF[jh] * y * get(node_share, nodes[jn], 0.0))


    # # nodal dispatch <= Available capacity in node
    # mod.ext[:constraints][:nodal_cap_limit] =  @constraint(mod, [jh in JH, jn in JN],
    #         g_nodal[jh,jn] ≤ y_nodal[jh,jn])
    
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH], 
        g[jh] <= y * AF[jh]) # capacity limit

    # offered capacity in capacity markets
    mod.ext[:constraints][:CM] = @constraint(mod, [jz=JZ], cap_cm[jz] <= σ_CM * PM[m][zones[jz]] * y)
    mod.ext[:constraints][:CM_sum] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) <= y)


    
    mod.ext[:constraints][:capcm_share] = @constraint(mod, [jn in JN], cap_cm_nodal[jn] == cap_cm[zone_of_idx[jn]]  * get(node_share, nodes[jn], 0.0))

    return mod

end