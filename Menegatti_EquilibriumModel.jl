###############################################################################
# Cross-Border participation: A false hope for fixing capacity market externalities?
# This code is currently written as a game theoretic equilibrium model but should be solved as a mixed complementarity problem (MCP).
###############################################################################
using JuMP
using Gurobi

model = Model(Gurobi.Optimizer)

###############################################################################
# 1. Parameters
###############################################################################

# -------------------------
# Sets
# -------------------------
zones = [:A, :B]                       # Electricity market zones
technologies = [:Base, :Mid, :Peak]    # Conventional generation technologies
renewables = [:Wind, :Solar]           # Renewable generation technologies
representative_days = 1:10
time_steps = 1:24

representative_day_weights = [1, 9, 9, 82, 11, 11, 180, 22, 22, 18]  # Example day weights

# -------------------------
# Cost/Price‐Related
# -------------------------
# Value of Lost Load (VOLL)
VOLL = Dict(:A => 3000, :B => 3000)  # €/MWh

# Variable costs for technologies (v_{c,k})
variable_costs = Dict(:Base => 36, :Mid => 53, :Peak => 76)  # €/MWh

# Annualized investment costs (I_{c,k})
investment_costs = Dict(
    :A => Dict(:Base => 180.0 + 1E4, :Mid => 100.0, :Peak => 70.0 + 1E4),
    :B => Dict(:Base => 180.0, :Mid => 100.0 + 1E4, :Peak => 70.0)
)

# We assume the user provides or calculates the *energy market price* λ_{c,t,day}
# For illustration, we define a placeholder structure:

# energy_price = dual of market clearing constraint

# -------------------------
# Capacity Market Prices
# -------------------------
# Local capacity‐market price (λ_c^{MC})
capacity_market_price = Dict(:A => 0.0, :B => 0.0)
# Foreign capacity‐market price (λ_c^{MC'}) 
# => e.g., zone A can also bid into B's capacity market and vice versa
capacity_market_price_foreign = Dict(:A => 0.0, :B => 0.0)

# You can set these to nonzero for whichever zone(s) have an active capacity market
# e.g., capacity_market_price[:A] = 100.0, capacity_market_price_foreign[:A] = 100.0, etc.

# -------------------------
# System Constraints
# -------------------------
capacity_margin = 0.10                # 10%
interconnection_capacity = 6000       # MW
maximum_capacity_price = 1_000_000    # €/MW

renewable_availability = Dict(
    :A => Dict(:Wind => 0.25, :Solar => 0.15),
    :B => Dict(:Wind => 0.30, :Solar => 0.20)
)

##############################################################################
# Energy Demand and Unserved Energy (Section 2.2.1)
##############################################################################
@variable(model, 0 <= unserved_energy[zone in zones, t in time_steps, day in representative_days],
          basename = "ENS_{c,t,day}") # Unserved Energy (MWh)

@variable(model, 0 <= energy_demand[zone in zones, t in time_steps, day in representative_days],
          basename = "d_{c,t,day}") # Energy Demand (MWh)

# We assume you have defined or loaded:
#    reference_demand[zone, t, day]  -> D_{c,t} (maximum or "reference" demand).

# ----------------------------------------------------------------------------
#  Objective (Consumer Surplus)
# ----------------------------------------------------------------------------
# We maximize the total consumer surplus, which is the sum over all zones,
# time steps, and representative days of:
#       (VOLL_c - λ_{c,t}) * d_{c,t}
# matching Section 2.2.1's approach, where λ_{c,t} is the electricity price.
# ----------------------------------------------------------------------------

@objective(model, Max, 
    sum( (VOLL[zone] - electricity_market_price[zone, t, day]) * energy_demand[zone, t, day]
         for zone in zones
         for t in time_steps
         for day in representative_days)
)

# ----------------------------------------------------------------------------
#  Constraints for Energy Demand (Equation (d.1), Section 2.2.1)
# ----------------------------------------------------------------------------
# Equation (d.1) states: 0 <= d_{c,t} <= D_{c,t},
# i.e., the served demand cannot exceed the reference demand, nor be negative.
@constraint(model, [zone in zones, day in representative_days, t in time_steps],
    0 <= energy_demand[zone, t, day] <= reference_demand[zone, t, day]
)

# In addition, unserved_energy_{c,t} = D_{c,t} - d_{c,t}.
@constraint(model, [zone in zones, day in representative_days, t in time_steps],
    unserved_energy[zone, t, day] == reference_demand[zone, t, day] - energy_demand[zone, t, day]
)


# -------------------------
# Generation & Capacity - Section 2.2.2
# -------------------------
@variable(model, 0 <= energy_generation[zone in zones, tech in technologies, t in time_steps, day in representative_days], 
          basename = "g_{c,k,t,day}", upperbound = 1.0e9)  # Energy Generation (MWh) - big number if needed for upper bound
@variable(model, 0 <= installed_capacity[zone in zones, tech in technologies],
          basename = "cP_{c,k}") # Installed Capacity (MW)

###############################################################################
# 3. Objective Function (Section 2.2.2, eqn. (g.1)–(g.2b))
###############################################################################
# In the text: 
#   max ∑_t [ (λ_{c,t} - v_{c,k}) * g_{c,k,t} ]
#       + λ_c^{MC} * cP_{c,k}^{MC}
#       + λ_c^{MC'} * cP_{c,k}^{MC'}
#       - I_{c,k} * cP_{c,k}
#
# Typically, this is summed over zones, technologies, time. We also incorporate
# representative_day_weights if we want to weight each day to reflect annual totals.

@objective(model, Max,
    sum(
        representative_day_weights[day] * (
            # Summation over hours: (energy price - variable cost) * generation
            sum( (energy_price[zone][day][t] - variable_costs[tech])
                 * energy_generation[zone, tech, t, day]
                 for t in time_steps )
            # Local capacity revenue
            + capacity_market_price[zone] * capacity_market_procurement[zone, tech]
            # Foreign capacity revenue
            + capacity_market_price_foreign[zone] *
              capacity_market_procurement_foreign[zone, tech]
            # Subtract investment cost
            - investment_costs[zone][tech] * installed_capacity[zone, tech]
        )
        for zone in zones
        for tech in technologies
        for day in representative_days
    )
)

###############################################################################
# 4. Constraints
###############################################################################
# (g.1) Generation bounded by installed capacity:
#       0 <= g_{c,k,t} <= cP_{c,k}
@constraint(model, [zone in zones, tech in technologies, t in time_steps, day in representative_days],
    energy_generation[zone, tech, t, day] <= installed_capacity[zone, tech]
)

# (g.2a) Local capacity‐market participation limited by installed capacity:
#        0 <= cP_{c,k}^{MC} <= cP_{c,k}
@constraint(model, [zone in zones, tech in technologies],
    capacity_market_procurement[zone, tech] <= installed_capacity[zone, tech]
)

# (g.2b) Foreign capacity‐market participation also limited by installed capacity:
#        0 <= cP_{c,k}^{MC'} <= cP_{c,k}
@constraint(model, [zone in zones, tech in technologies],
    capacity_market_procurement_foreign[zone, tech] <= installed_capacity[zone, tech]
)

###############################################################################
# Section 2.2.3: Renewable Generators
# Each renewable generator j ∈ {Wind, Solar} in zone c decides on g_{c,j,t}.
###############################################################################
@variable(model, 0 <= renewable_generation[zone in zones, ren in renewables,
          day in representative_days, t in time_steps], basename = "g_{c,j,d,t}") # Renewable Generation (MWh)
@variable(model, 0 <= installed_capacity_renewable[zone in zones, ren in renewables],
          basename = "cP_{c,ren}") # Installed Renewable Capacity (MW)

# -------------------------
# 1. Availability factors
# -------------------------

@expression(model, renewable_availability_factor[zone in zones, ren in renewables,
                                             day in representative_days, t in time_steps],
    renewable_availability[zone][ren] * installed_capacity_renewable[zone, ren]
)

# 2. Objective: Maximize total revenue from the energy market
#    According to Section 2.2.3, eqn. (no explicit label in your text, 
#    but the formula is: max ∑_t [ λ_{c,t} * g_{c,j,t} ] )
#    We multiply each day’s contribution by `representative_day_weights[day]`
#    so that the solution reflects the annual (or seasonal) total.
@objective(model, Max, renewable_revenue, 
    sum(
        representative_day_weights[day] *
        energy_price[zone][day][t] * renewable_generation[zone, ren, day, t]
        for zone in zones
        for ren in renewables
        for day in representative_days
        for t in time_steps
    )
)

# 3. Constraints for Renewable Generators
#    Equation in your text: 0 <= g_{c,j,t} <= AF_{c,j,t} * cP_{c,j}
#    AF_{c,j,t} is the availability factor, cP_{c,j} is the installed capacity.
@constraint(model, [zone in zones, ren in renewables, day in representative_days, t in time_steps],
    renewable_generation[zone, ren, day, t] 
    <= renewable_availability_factor[zone, ren, day, t]
)

###############################################################################
# 2.2.4 Capacity market demand
###############################################################################

# A reference maximum capacity price in each bidding zone:
#   Λ_MAX^(CM,c),   Λ_MAX^(CM,c')
@parameter Λ_max_local = 1.0e6     # for zone c
@parameter Λ_max_foreign = 1.0e6   # for zone c'

# The actual capacity-market clearing prices:
#   λ_c^(CM),   λ_c'^(CM,c)
@parameter λ_c_local   = 120.0  # example
@parameter λ_c_foreign = 110.0  # example

# The capacity demand targets appearing in constraints (d.3) and (d.4).
#   D_c'^(CM,c) and D_(c+c')^(CM,c')
# In your text, D_c'^(CM,c) is the foreign capacity target from zone c',
# while D_(c+c')^(CM,c') is the total capacity target that includes
# cross-border availability.
@parameter D_cprime_local = 5000.0   # = D_c'^(CM,c)  (limit on foreign demand)
@parameter D_cpluscprime_foreign = 8000.0  # = D_(c+c')^(CM,c') (limit on total)

###############################################################################
# DECISION VARIABLES
###############################################################################
# d_c^(CM,c)   : capacity demanded by zone c in its own capacity market
# d_c^(CM,c')  : capacity demanded by zone c in the foreign zone c' capacity market
@variable(model, 0 <= d_local,    basename="d_c^(CM,c)")
@variable(model, 0 <= d_foreign,  basename="d_c^(CM,c')")
###############################################################################
# OBJECTIVE
#
#   max  (Λ_MAX^(CM,c)  -  λ_c^(CM))   * d_c^(CM,c)
#      + (Λ_MAX^(CM,c') -  λ_c'^(CM,c)) * d_c^(CM,c')
###############################################################################
@objective(model, Max, 
    (Λ_max_local   - λ_c_local)   * d_local   +
    (Λ_max_foreign - λ_c_foreign) * d_foreign
)

###############################################################################
# CONSTRAINTS
###############################################################################
# (d.3):  0 <= d_c^(CM,c') <= D_c'^(CM,c)
@constraint(model, constraint_d3, d_foreign <= D_cprime_local)

# (d.4):  0 <= d_c^(CM,c) + d_c^(CM,c') <= D_(c+c')^(CM,c')
@constraint(model, constraint_d4, d_local + d_foreign <= D_cpluscprime_foreign)
###############################################################################


###############################################################################
# Decision Problem for Interconnector
###############################################################################
# Objective for Interconnector: Maximize Congestion Revenue
@variable(model, 0 <= flow_c_to_cprime[t in time_steps] <= interconnection_capacity)
@variable(model, 0 <= flow_cprime_to_c[t in time_steps] <= interconnection_capacity)

# Possibly you enforce one direction at a time, or allow net flows, etc.

@objective(model, Max,
    sum((lambda_cprime[t] - lambda_c[t]) * flow_c_to_cprime[t]
      + (lambda_c[t] - lambda_cprime[t]) * flow_cprime_to_c[t]
        for t in time_steps)
)


###############################################################################
# Section 2.2.6 - Market Clearing
###############################################################################

# Suppose you already have:
#   - energy_generation[c, k, t, day]  for conventional gen (MW or MWh at hour t)
#   - renewable_generation[c, j, t, day] for renewable gen (MW or MWh at hour t)
#   - cross_border_flow[from_zone, to_zone, t, day] for energy flows (MWh)
#   - energy_demand[c, t, day] for load (MWh)
#   - unserved_energy[c, t, day] (optional)

# ---------------  (c.1) Energy Market Clearing  ---------------
# ∑_k g_{c,k,t} + ∑_j g_{c,j,t} + f_{c'->c,t} - f_{c->c',t} = d_{c,t}
#
# If you allow unserved demand, then the equation might read:
# ∑_k g_{c,k,t} + ∑_j g_{c,j,t} + f_{c'->c,t} - f_{c->c',t} + unserved_{c,t} = d_{c,t}
# which is an extension of (c.1).

@constraint(model, [zone in zones, t in time_steps, day in representative_days],
    sum(energy_generation[zone, tech, t, day] for tech in technologies)
  + sum(renewable_generation[zone, ren, t, day] for ren in renewables)
  + sum(cross_border_flow[other_zone, zone, t, day] for other_zone in zones if other_zone != zone)
  - sum(cross_border_flow[zone, other_zone, t, day] for other_zone in zones if other_zone != zone)
  + unserved_energy[zone, t, day]        # <= only if you want it
  == energy_demand[zone, t, day]
)

###############################################################################
# Section 2.2.6 - Capacity Market Clearing
###############################################################################
#
# (c.2)  ∑_k [ cP_{c,k}^{CM,c} ] = d_{c}^{CM,c}
# (c.3)  ∑_k [ cP_{c,k}^{CM,c'} ] = d_{c'}^{CM,c}

# Decision variables for capacity procurement:
#   capacity_market_procurement[c, k]       => cP_{c,k}^{CM,c}
#   capacity_market_procurement_foreign[c, k] => cP_{c,k}^{CM,c'}
# Demand parameters:
#   capacity_market_demand[c]               => d_{c}^{CM,c}
#   foreign_capacity_demand[c]              => d_{c'}^{CM,c} (the foreign demand c is fulfilling in zone c')

# ---------------  (c.2) Local capacity clearing ---------------
@constraint(model, [zone in zones],  # eqn. (c.2)
    sum(capacity_market_procurement[zone, tech] for tech in technologies)
    == capacity_market_demand[zone]
)

# ---------------  (c.3) Foreign capacity clearing ---------------
@constraint(model, [zone in zones],  # eqn. (c.3)
    sum(capacity_market_procurement_foreign[zone, tech] for tech in technologies)
    == foreign_capacity_demand[zone]
)

# That is the direct statement that the total capacity c buys from foreign zone c'
# equals the foreign capacity demand.

###############################################################################


for t in time_steps
    for c in zones
        for c_prime in zones
            if c != c_prime
                # Decrease f_{c->c',t} if c has unserved energy
                Δf_c_cprime = -min(ens[c, t], f[c, c_prime, t])
                f[c, c_prime, t] += Δf_c_cprime
                ens[c, t]       += Δf_c_cprime   # ens decreases if Δf is negative

                # Decrease f_{c'->c,t} if c' has unserved energy
                Δf_cprime_c = -min(ens[c_prime, t], f[c_prime, c, t])
                f[c_prime, c, t] += Δf_cprime_c
                ens[c_prime, t]  += Δf_cprime_c

                # Then recalc or adjust demand if needed 
                # (or keep track if your model tracks d_c,t changes)
            end
        end
    end
end


