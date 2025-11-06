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
    λ_EOM = mod.ext[:parameters][:λ_EOM]            # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]            # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]            # rho-value in ADMM related to EOM auctions
    cap_bar = mod.ext[:parameters][:cap_bar]        # element in ADMM penalty term related to capacity markets
    ρ_CM = mod.ext[:parameters][:ρ_CM]              # rho-value in ADMM related to capacity market auctions
    λ_CM = mod.ext[:parameters][:λ_CM]              # CM prices
    σ_CM = mod.ext[:parameters][:σ_CM]              # 1 if capacity markets are active, 0 otherwise
    PM = mod.ext[:parameters][:participation_matrix]                # participation matrix
    node_share = mod.ext[:parameters][:node_share_vec]              # Vector{Float64}


    # Investment parameters
    I = mod.ext[:parameters][:I]                    # investment cost
    max_cap = mod.ext[:parameters][:max_cap]        # maximum installable capacity of the generator
    y_init = mod.ext[:parameters][:C]               # existing capacity
    A = mod.ext[:parameters][:A]                    # cost function parameter - quadratic cost
    B = mod.ext[:parameters][:B]                    # cost function parameter - marginal cost

    # Variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], lower_bound=0, base_name="generation")                        # generation
    y = mod.ext[:variables][:y] = @variable(mod, lower_bound=0, base_name="InstalledCapacity")                     # installed capacity
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], lower_bound = 0, base_name = "capacity_offered")    # capacity offered in capacity markets

    # Objective => minimize GenCo costs
    mod.ext[:objective] = @objective(mod, Min,
        + sum(W[jh] * A/2*g[jh]^2 for jh in JH)                                             # cost function for generation
        + sum(W[jh] * B*g[jh] for jh in JH)
        - sum(W[jh] * λ_EOM[jh]*g[jh] for jh in JH)                                         # weighted revenue from EOM
        + I * y                                                                             # annualized investment cost
        - σ_CM * sum(λ_CM[home_zone] * cap_cm[jz] for jz in JZ)                            # capacity market revenue - uniform price auction
        + sum(W[jh] * ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)                           # ADMM penalty term for EOM clearing with weights
        + σ_CM * ρ_CM[home_zone]/2 * (sum(cap_cm[jz] for jz in JZ) - cap_bar[home_zone])^2  # ADMM penalty term for CM
    )

    # Nodal available capacity expression
    mod.ext[:expressions][:y_nodal] = @expression(mod, [jh in JH, jn in JN],
        y_init * node_share[jn] * AF[jh]
    )

    # generation ≤ available capacity
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH], 
        g[jh] <= (y + y_init) * AF[jh])

    # Renewable capacity constraint
    mod.ext[:constraints][:ren_cap] = @constraint(mod, 
        y + y_init <= max_cap)

    # Capacity market constraints
    mod.ext[:constraints][:CM] = @constraint(mod, [jz=JZ], 
        cap_cm[jz] <= σ_CM * PM[m][zones[jz]] * (y + y_init)) # capacity offered cannot exceed installed capacity * participation factor [0,1]

    mod.ext[:constraints][:CM_sum] = @constraint(mod, 
        sum(cap_cm[jz] for jz in JZ) <= (y + y_init))

    return mod
end