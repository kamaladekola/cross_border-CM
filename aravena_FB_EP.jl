using JuMP, Gurobi
using Plots
# ----------------------- DATA SECTION ---------------------------------
const ZONES  = [:A, :B, :C]
const NODES  = [:n1, :n2, :n3, :n4]
const zone_of = Dict(:n1=>:A, :n2=>:A, :n3=>:B, :n4=>:C)

# branches  (from, to, thermal limit F [MW])
const BRANCHES = [
    (:n1, :n2, 100.0),   # ℓ₁₂   (A–A)
    (:n2, :n3, 2000.0),   # ℓ₂₄   (A–B)
    (:n3, :n4, 2000.0),   # ℓ₄₃   (B–C)
    (:n4, :n1, 2000.0)    # ℓ₃₁   (C–A)
]

# the 4×4 PTDF matrix (rows = branches 1–4, cols = nodes n1–n4)
const PTDF = [
     0.5  -0.25   0.25   0.0;   # l1: 1–2
     0.5   0.75   0.25   0.0;   # l2: 2–4
    -0.5  -0.25  -0.75   0.0;   # l3: 4–3
    -0.5  -0.25   0.25   0.0    # l4: 3–1
]

# generator capacities
const GENERATORS = [
    (:n1, 500.0),
    (:n2,  200.0),
    (:n3,  300.0),
    (:n4,  500.0)
]

const COST = [5.0, 70.0, 15.0, 40.0]    # €/MWh for gens 1–4

# fixed demands
const DEMAND = Dict(:n1=>0.0, :n2=>300.0, :n3=>0.0, :n4=>300.0)
# ---------------------------------------------------------------------

function feasible_net_position(p_vec::Vector{<:Real}; atol=1e-6)
    @assert length(p_vec) == length(ZONES)
    p_val = Dict(ZONES .=> p_vec)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    # 1) dispatch fractions v_g ∈ [0,1]
    @variable(model, 0 <= v[1:length(GENERATORS)] <= 1)

    # 2) branch flows f_ℓ
    @variable(model, f[1:length(BRANCHES)])

    # 3) zonal balance
    for z in ZONES
        @constraint(model,
            sum(GENERATORS[i][2]*v[i]
                for i in 1:length(GENERATORS)
                if zone_of[GENERATORS[i][1]] == z)
            - p_val[z]
            == sum(DEMAND[n] for n in NODES if zone_of[n]==z)
        )
    end

    # 4) relate f_ℓ to nodal net injections via PTDF
    #    net injection at node n = Σ_{g@node n} Q_g·v_g  - DEMAND[n]
    for ℓ in eachindex(BRANCHES)
        @constraint(model,
            f[ℓ] == sum(
                PTDF[ℓ, j] * (
                    sum(
                        GENERATORS[i][2]*v[i]
                        for i in eachindex(GENERATORS)
                        if GENERATORS[i][1] == NODES[j]
                    )
                  - DEMAND[NODES[j]]
                )
                for j in eachindex(NODES)
            )
        )
        # 5) thermal limits
        _,_,F = BRANCHES[ℓ]
        @constraint(model, -F <= f[ℓ] <= F)
    end

    @objective(model, Min, sum( COST[i] * GENERATORS[i][2] * v[i] 
    for i in eachindex(GENERATORS)))

    optimize!(model)

    return termination_status(model)==MOI.OPTIMAL ||
           primal_status(model)==MOI.FEASIBLE_POINT
end

# --------------------------- Feasibility test -------------------------------------
    
for pair in ([-50, 100, -50],
    [0, 100, -100],
    [180, -90, -90] ,
)  
    p = pair
    ok = feasible_net_position(p)
    println(" p=",pair,"   → ", ok ? "feasible" : "INFEASIBLE")
end