function solve_generator_agent!(mod::Model, m::String, zones::Vector{String})

    zone, _ = parse_agent_name(m)
    home_zone = findfirst(isequal(zone), zones)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    W = mod.ext[:parameters][:w]

    # Extract parameters
    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM]        # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]        # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]        # rho-value in ADMM related to EOM auctions
    cap_bar = mod.ext[:parameters][:cap_bar]    # element in ADMM penalty term related to capacity markets
    ρ_CM = mod.ext[:parameters][:ρ_CM]          # rho-value in ADMM related to capacity markets
    λ_CM = mod.ext[:parameters][:λ_CM]          # CM prices
    I = mod.ext[:parameters][:I]                # investment cost
    y_init   = mod.ext[:parameters][:C]         # existing capacity
    σ_CM = mod.ext[:parameters][:σ_CM]          # 1 if capacity markets are active, 0 otherwise


    # Create variables
    g = mod.ext[:variables][:g]  
    y = mod.ext[:variables][:y]
    cap_cm = mod.ext[:variables][:cap_cm]


    # Objective => minimize GenCo costs
    mod.ext[:objective] = @objective(mod, Min,
        + sum(W[jh] * A/2*g[jh]^2 for jh in JH)                                             # cost function for generation
        + sum(W[jh] * B*g[jh] for jh in JH)
        - sum(W[jh] * λ_EOM[jh]*g[jh] for jh in JH)                                         # weighted revenue from EOM
        + I * y                                                                             # annualized investment cost
        - σ_CM * sum(λ_CM[home_zone] * cap_cm[jz] for jz in JZ)                             # capacity market revenue - uniform price auction
        + sum(W[jh] * ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)                           # ADMM penalty term for EOM clearing with weights
        + σ_CM * ρ_CM[home_zone]/2 * (sum(cap_cm[jz] for jz in JZ) - cap_bar[home_zone])^2  # ADMM penalty term for CM
    )


    optimize!(mod);

    return mod
end

