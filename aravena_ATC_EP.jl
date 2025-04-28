"""
ATC-EP: Available-Transfer-Capacity with Exact Projection
Tests if a net position profile p = (p_A, p_B, p_C) belongs to the ATC-EP
"""
# (:A, :B), (:B, :C),(:C, :A) 
ATC = [[170.0, -510.0], [340.0, -340.0], [510.0, -170.0]] # ATC limits for each interconnector
# ----------------------- IMPORTS --------------------------------------
using JuMP, Gurobi

# ----------------------- DATA SECTION ---------------------------------
const ZONES  = [:A, :B, :C]
const NODES  = [:n1, :n2, :n3, :n4]
const zone_of = Dict(:n1=>:A, :n2=>:A, :n3=>:B, :n4=>:C)

# generator capacities
const GENERATORS = [
    (:n1, 500.0),  # gen at node n1 (zone A)
    (:n2, 200.0),  # n2 (A)
    (:n3, 300.0),  # n3 (B)
    (:n4, 500.0)   # n4 (C)
]

const COST = [5.0, 70.0, 15.0, 40.0]    # €/MWh for gens 1–4


# fixed demand at each node
const DEMAND = Dict(:n1=>0.0, :n2=>300.0, :n3=>0.0, :n4=>300.0)

# branches = (from, to, thermal limit F [MW])
const BRANCHES = [
    (:n1, :n2, 100.0),   # ℓ₁₂   (A–A) 
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

# Interconnectors T ⊆ Z×Z (only neighboring zones)
const TCONNECT = [(:A,:B), (:B,:C), (:C,:A)]




# STEP 2: Feasibility check per Eq 13
function atc_ep_feasible(p_vec::Vector{<:Real})
    @assert length(p_vec) == length(ZONES)
    p_val = Dict(ZONES .=> p_vec)

    atc_plus = Dict()
    atc_minus = Dict()
    for (i, t) in enumerate(TCONNECT)
        atc_plus[t] = ATC[i][1]       # First value is positive capacity
        atc_minus[t] = -ATC[i][2]     # Second value is negative capacity
    end

    m = Model(Gurobi.Optimizer)
    set_silent(m)

    # 1) dispatch fractions
    @variable(m, 0 <= v[i=1:length(GENERATORS)] <= 1)

    # 2) *all* branch flows
    @variable(m, f[ℓ=1:length(BRANCHES)])

    # 3) cross‐border flows as a separate var, tied to f
    @variable(m, e[t in TCONNECT])
    for (ℓ,(i,j,F)) in enumerate(BRANCHES)
        # if (i,j) is A–B then connect f[ℓ] == e[(:A,:B)], etc.
        if (zone_of[i],zone_of[j]) in TCONNECT
            @constraint(m, f[ℓ] ==  e[(zone_of[i],zone_of[j])])
        elseif (zone_of[j],zone_of[i]) in TCONNECT
            @constraint(m, f[ℓ] == -e[(zone_of[j],zone_of[i])])
        end
        # thermal limits on every line
        @constraint(m, -F <= f[ℓ] <= F)
    end

    # 4) zonal balance and PTDF for *every* branch
    for z in ZONES
        @constraint(m,
            sum(GENERATORS[i][2]*v[i]
                for i in 1:length(GENERATORS)
                if zone_of[GENERATORS[i][1]] == z)
            - p_val[z]
            == sum(DEMAND[n] for n in NODES if zone_of[n]==z)
        )
    end
    for ℓ in eachindex(BRANCHES)
        @constraint(m,
            f[ℓ] == sum(
                PTDF[ℓ,j] * (
                    sum(GENERATORS[i][2]*v[i]
                        for i in eachindex(GENERATORS)
                        if GENERATORS[i][1] == NODES[j])
                  - DEMAND[NODES[j]]
                )
                for j in eachindex(NODES)
            )
        )
    end

    # 5) ATC bounds on the e[t] variables
    for t in TCONNECT
        @constraint(m, -atc_minus[t] <= e[t] <= atc_plus[t])
    end

    
    @objective(m, Min, sum( COST[i] * GENERATORS[i][2] * v[i] 
    for i in eachindex(GENERATORS)))

    optimize!(m)
    return termination_status(m) == MOI.OPTIMAL ||
           primal_status(m)    == MOI.FEASIBLE_POINT
end



println("\n=== Feasibility tests (Eq 13) ===")

for pair in ([150, -100, -50],
             [0, 100, -100],
             [180, -90, -90] ,
    )  
p = pair
ok = atc_ep_feasible(p)
println(" p=",pair,"   → ", ok ? "feasible" : "INFEASIBLE")
end