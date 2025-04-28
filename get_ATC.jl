"""
Implements the ATC-EP method from Aravena et al. (2021) "Transmission Capacity 
Allocation in Zonal Electricity Markets".(Eqs 9-12)
"""


# using Pkg
# Pkg.add("Ipopt") 
using JuMP, Gurobi, Ipopt 

const ZONES  = [:A, :B, :C]
const NODES  = [:n1, :n2, :n3, :n4]
const zone_of = Dict(:n1=>:A, :n2=>:A, :n3=>:B, :n4=>:C)

# generator capacities
const GENERATORS = [
    (:n1, 1000000000.0),  # gen at node n1 (zone A)
    (:n2, 1000000000.0),  # n2 (A)
    (:n3, 1000000000.0),  # n3 (B)
    (:n4, 1000000000.0)   # n4 (C)
]


# const COST = [5.0, 70.0, 15.0, 40.0]    # €/MWh for gens 1–4


# fixed demand at each node
const DEMAND = Dict(:n1=>500.0, :n2=>500.0, :n3=>500.0, :n4=>500.0)

#  thermal limit F [MW]
const BRANCHES = [
    (:n1, :n2, 500.0),   # ℓ₁₂   (A–A)
    (:n2, :n3, 2000.0),   # ℓ₂₄   (A–B)
    (:n3, :n4, 2000.0),   # ℓ₄₃   (B–C)
    (:n4, :n1, 2000.0)    # ℓ₃₁   (C–A)
]


# node‐to‐line PTDF[ℓ,n]
const PTDF = [
     0.5  -0.25   0.25   0.0;   # l1: 1–2
     0.5   0.75   0.25   0.0;   # l2: 2–4
    -0.5  -0.25  -0.75   0.0;   # l3: 4–3
    -0.5  -0.25   0.25   0.0    # l4: 3–1
]

# Interconnectors
const TCONNECT = [(:A,:B), (:B,:C), (:C,:A)]

# Create  optimization model
# atc_model = Model(Gurobi.Optimizer)
atc_model = Model(Ipopt.Optimizer)
set_silent(atc_model)

# ATC variables
@variable(atc_model, 0 ≤ atc_plus[t in TCONNECT])
@variable(atc_model, 0 ≤ atc_minus[t in TCONNECT])

# Eq.(10)
@constraint(atc_model, [t in TCONNECT], -atc_minus[t] <= atc_plus[t])


# 2^3 vertices:
num_interconnects = length(TCONNECT)
signatures = collect(Iterators.product(fill([-1,1], num_interconnects)...))
nverts     = length(signatures)

# variables for each rectangle-vertex ν
@variable(atc_model, 0 <= v[ν=1:nverts, i=1:length(GENERATORS)] <= 1) # utilisation factor of generator i when the vertex ν exchanges the corresponding ATC.
@variable(atc_model,     f[ν=1:nverts, ℓ=1:length(BRANCHES)]) # realised flow (MW) on physical line ℓ at vertex ν.
@variable(atc_model,     pz[ν=1:nverts, z in ZONES]) # zonal net position (MW) at vertex ν

# For each νertex, enforce DC-flow feasibility
for (ν,sig) in enumerate(signatures)

    # e = exchange on each interconnector
    e = Dict{Tuple{Symbol,Symbol},JuMP.AffExpr}()
    for (k,t) in enumerate(TCONNECT)
        e[t] = sig[k] == +1 ?  atc_plus[t] :  -atc_minus[t]
    end



    # zonal Demand-Supply balance
    for z in ZONES
        @constraint(atc_model,
          sum(GENERATORS[i][2]*v[ν,i] for i in eachindex(GENERATORS) if zone_of[GENERATORS[i][1]]==z) - pz[ν,z]
          == sum(DEMAND[n] for n in NODES if zone_of[n]==z)
        )

        @constraint(atc_model,
        pz[ν,z] == sum(e[t] for t in TCONNECT if t[1]==z) - sum(e[t] for t in TCONNECT if t[2]==z)) # net position for zone z

        # @constraint(atc_model, sum(pz[ν,z] for z in ZONES) == 0) # net position balance
    end

    
    # Nodal balance
    for ℓ in eachindex(BRANCHES)
        (i, j, F) = BRANCHES[ℓ]
        # Calculate flows using PTDF
        @constraint(atc_model,
            f[ν,ℓ] == sum( PTDF[ℓ,j] * (sum(GENERATORS[g][2]*v[ν,g] for g in eachindex(GENERATORS) if GENERATORS[g][1] == NODES[j])
                  - DEMAND[NODES[j]]) for j in eachindex(NODES))
        )
        @constraint(atc_model, -F <= f[ν,ℓ] <= F) # Thermal limits
    end
end


# @objective(atc_model, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT))
@NLobjective(atc_model, Max, prod(atc_plus[t] + atc_minus[t] for t in TCONNECT))


optimize!(atc_model)




ATCplus  = value.(atc_plus)
ATCminus = value.(atc_minus)

println("=== ATC-EP results (Eqs 9–12) ===")
for t in TCONNECT
  println("  ", t,
          "   +", round(ATCplus[t], digits=1),
          "   -", round(ATCminus[t], digits=1))
end