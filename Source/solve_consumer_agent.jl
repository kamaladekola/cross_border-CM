function solve_consumer_agent!(mod::Model)
   # Extract sets
   JH = mod.ext[:sets][:JH]

   # Extract parameters
   λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
   g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
   ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

   # Create variables
   g = mod.ext[:variables][:g]  
  
   # Objective 
   mod.ext[:objective] = @objective(mod, Min,
       - sum(λ_EOM[jh]*g[jh] for jh in JH)
       + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)
   )

    optimize!(mod);

    return mod
end