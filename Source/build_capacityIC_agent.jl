function build_capacityIC_agent!(mod::Model)
    JZ = mod.ext[:sets][:JZ]
    JL = mod.ext[:sets][:JL]
    # JS = mod.ext[:sets][:JS]
    # JB = mod.ext[:sets][:JB]


    PTDF = mod.ext[:parameters][:PTDF]          # (nLines × nZones)
    RAM = mod.ext[:parameters][:RAM]            # (nLines)
    coupling = mod.ext[:parameters][:coupling]      # "FB" or "ATC"

    # border_names = mod.ext[:sets][:border_names] 
    # MEC  = mod.ext[:parameters][:MEC]           # (nZones)
    # from_zone  = mod.ext[:parameters][:from_zone]
    # to_zone    = mod.ext[:parameters][:to_zone]

    # ADMM penalty parameters for capacity market
    cap_bar = mod.ext[:parameters][:cap_bar] # ADMM consensus variables
    λ_CM    = mod.ext[:parameters][:λ_CM]
    ρ_CM    = mod.ext[:parameters][:ρ_CM]

    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], base_name = "zonal_capacity_netpos")


    # Objective
    mod.ext[:objective] =  @objective(mod, Min,
    - sum(λ_CM[jz] * cap_cm[jz] for jz in JZ) 
    + sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
)

    # Constraints
    mod.ext[:constraints][:CapacityNP] = @constraint(mod, sum(cap_cm[jz] for jz in JZ) == 0)

    if coupling == "FB"

        #   mod.ext[:constraints][:FBMC] = @constraint(mod,[jh=JH, jl=JL], -RAM[jl] <= sum(PTDF[jl, jz] * g[jh, jz] for jz in JZ) <= RAM[jl]) # scenario based dispatch to simulate different scarcity moment in different zones

        mod.ext[:constraints][:CapacityIC] = @constraint(mod, [jl in JL], -RAM[jl] <= sum(PTDF[jl, jz] * cap_cm[jz] for jz in JZ ) <= RAM[jl]) # flow based market coupling for capacity offered in capacity markets
        # include clearing -> supply = demand (to ensure feasible dispatch)
        # g[jh, jz] <= cap_cm  # dispatched electricity across interconnector must be less than forigen capacity offer.

    
    elseif coupling == "ATC" 
    # each border is independently constrained by MEC 
    # exact projection onto the domain of FBMC to ensure network feasibility
     #   mod.ext[:constraints][:FBMC] = @constraint(mod,[jh=JH, jl=JL], -RAM[jl] <= sum(PTDF[jl, jz] * g[jh, jz] for jz in JZ) <= RAM[jl]) # scenario based dispatch to simulate different scarcity moment in different zones
     # include clearing -> supply = demand (to ensure feasible dispatch)
     # g[jh, jz] <= cap_cm  # dispatched electricity across interconnector must be less than forigen capacity offer.
    end
    return mod
end


