function build_interconnector_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract parameters
    λ1 = mod.ext[:parameters][:λ1]
    λ2 = mod.ext[:parameters][:λ2]
    IC = mod.ext[:parameters][:IC]
    IC_r = mod.ext[:parameters][:IC_r]

    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    # Create variables
    # -IC < f_{z1->z2,t} >= IC  (net flow -> negative for export)
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], lower_bound=-IC_r, upper_bound=IC, base_name="generation") # net flow -> negative for export (i.e if power flows from c to c')

    # Objective: Maximize congestion revenue => sum( (λ2[t] - λ1[t]) * f[t])
    mod.ext[:objective] = @objective(mod, Min,
    - sum((λ2[jh] - λ1[jh]) * g[jh] for jh in JH)       # when g > 0, flow is from zone  1 to zone 2.
    + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
)


    return mod
end


# as defined in Menegatti et al.

# enforce one direction at a time, or allow net flows?

# extend to include flow based market coupling constraints?


