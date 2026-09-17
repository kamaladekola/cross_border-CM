

# import Pkg
# Pkg.add("MathOptInterface")
using JuMP
using Gurobi
using DataFrames, CSV, YAML
import MathOptInterface as MOI

# Helper functions
function get_any(d::Dict, keys::AbstractVector{<:AbstractString})
    for k in keys
        if haskey(d, k);            return d[k]            end
        if haskey(d, lowercase(k)); return d[lowercase(k)] end
        if haskey(d, uppercase(k)); return d[uppercase(k)] end
    end
    error("Missing keys $(join(keys, ", ")) in dictionary.")
end
get_any(d::Dict, key::AbstractString) = get_any(d, [key])

function safe_col(df::DataFrame, name::Union{String,Symbol})
    s = Symbol(name)
    @assert s ∈ propertynames(df) "Missing column $(s) in $(nameof(parentmodule(df))) DataFrame"
    return df[!, s]
end

struct Sets
H::Int
Z::Int
N::Int
I::Int
L::Int
S::Int
JH::UnitRange{Int}
JZ::UnitRange{Int}
JN::UnitRange{Int}
JI::UnitRange{Int}
JL::UnitRange{Int}
JS::UnitRange{Int}
zones::Vector{String}
nodes::Vector{Symbol}
techs::Vector{Symbol}
end

struct Maps
    zone_of_node::Vector{Int} # length N, zone index per node
    zone_to_idx::Dict{String,Int}
    tech_to_idx::Dict{Symbol,Int}
    zone_nodes_idx::Vector{Vector{Int}} # per z, list of node indices
end

struct Params
    W::Matrix{Float64} # weights per t
    D_zonal::Matrix{Float64} # H × Z (demand)
    AV::Array{Float64,3} # H × I × N (availability factor)
    y_node::Matrix{Float64} # I × N (existing capacity)
    MC::Matrix{Float64} # I × Z (marginal cost)
    max_cap::Matrix{Float64} # I × Z (maximum capacity)
    A::Matrix{Float64} # I × Z (availability)
    IC::Matrix{Float64} # I × Z (investment cost)
    S::Matrix{Float64} # Z × N (demand share per zone-node)
    PTDF::Matrix{Float64} # L × N (power transfer distribution factors)
    Fmax::Vector{Float64} # L (thermal limits)

    # CM parameters
    margins::Vector{Float64} # capacity margin per zone
    CD_ref::Vector{Float64} # capacity demand per zone
    scarcity_matrix::Matrix{Float64} # S × Z demand in scarcity scenarios
end

function read_inputs(base_dir::AbstractString)
    data = YAML.load_file(joinpath(base_dir, "Input", "config.yaml"))
    load = CSV.read(joinpath(base_dir, "Input", "load.csv"), DataFrame; delim = ";")
    pv = CSV.read(joinpath(base_dir, "Input", "pv.csv"), DataFrame; delim = ";")
    wind_on = CSV.read(joinpath(base_dir, "Input", "wind_onshore.csv"), DataFrame; delim = ";")
    nptdf = CSV.read(joinpath(base_dir, "Input", "nodal_ptdf.csv"), DataFrame; delim = ";")
    lines = CSV.read(joinpath(base_dir, "Input", "lines.csv"), DataFrame; delim = ";")
    scarcity_df = CSV.read(joinpath(base_dir, "Input", "scarcity.csv"), DataFrame; delim = ";")


    H = Int(data["General"]["nTimesteps"])
    weights = joinpath(base_dir, "Input", "weights.csv")
    W = if isfile(weights)
    df = CSV.read(weights, DataFrame; delim = ";")
    :weights ∈ propertynames(df) ? Float64.(df[!, :weights][1:H]) : ones(Float64, H)
    else
    ones(Float64, H)
    end


    return (data=data, load=load, pv=pv, wind_on=wind_on,
    nptdf=nptdf, lines=lines, weights=W, scarcity=scarcity_df)
end

function read_inputs(base_dir::AbstractString)
    data = YAML.load_file(joinpath(base_dir, "Input", "config.yaml"))
    load = CSV.read(joinpath(base_dir, "Input", "load.csv"), DataFrame; delim = ";")
    pv = CSV.read(joinpath(base_dir, "Input", "pv.csv"), DataFrame; delim = ";")
    wind_on = CSV.read(joinpath(base_dir, "Input", "wind_onshore.csv"), DataFrame; delim = ";")
    nptdf = CSV.read(joinpath(base_dir, "Input", "nodal_ptdf.csv"), DataFrame; delim = ";")
    lines = CSV.read(joinpath(base_dir, "Input", "lines.csv"), DataFrame; delim = ";")
    scarcity_df = CSV.read(joinpath(base_dir, "Input", "scarcity.csv"), DataFrame; delim = ";")

    H = Int(data["General"]["nTimesteps"])
    zones = collect_zones(data)
    Z = length(zones)
    
    weights = joinpath(base_dir, "Input", "weights.csv")
    W = if isfile(weights)
        df = CSV.read(weights, DataFrame; delim = ";")
        # Read zone-specific columns or use uniform weights
        W_matrix = zeros(Float64, H, Z)
        for (zidx, z) in enumerate(zones)
            col_sym = Symbol(z)  # Changed from Symbol("(z)") to Symbol(z)
            if col_sym ∈ propertynames(df)
                W_matrix[:, zidx] = Float64.(df[!, col_sym][1:H])
            # elseif :weights ∈ propertynames(df)
            #     # Fallback: use single weights column for all zones
            #     W_matrix[:, zidx] = Float64.(df[!, :weights][1:H])
            else
                error("No weight column found for zone $(z). Expected '$(z)' or 'weights' in weights.csv")
            end
        end
        W_matrix
    else
        error("weights.csv file not found")  # Also fixed the error message here
    end

    return (data=data, load=load, pv=pv, wind_on=wind_on,
    nptdf=nptdf, lines=lines, weights=W, scarcity=scarcity_df)
end

function collect_zones(data)
    return sort!(String.(collect(keys(data["Consumers"]))))
end

function collect_techs(data)

    tech_names = String[]
    for (_zone, gens) in data["Generators"]
        for k in keys(gens)
            push!(tech_names, String(k))
        end
    end
    tech_names = unique(tech_names)
    sort!(tech_names)
    return Symbol.(tech_names)
end

function build_sets_and_maps(data, lines::DataFrame, scarcity_df::DataFrame)
    zones = collect_zones(data)
    techs = collect_techs(data)
    zone_nodes = Dict{String,Vector{Symbol}}(z => Symbol.(data["Network"]["ZoneMap"][z]) for z in zones)
    nodes = unique(vcat([Symbol.(data["Network"]["ZoneMap"][z]) for z in zones]...))


    H = Int(data["General"]["nTimesteps"]) ; Z = length(zones) ; N = length(nodes)
    I = length(techs) ; L = nrow(lines)
    S = nrow(scarcity_df)

    zone_to_idx = Dict(z => i for (i, z) in enumerate(zones))
    tech_to_idx = Dict(techs[i] => i for i in 1:I)

    zone_of_node = similar(Vector{Int}(), N)
    resize!(zone_of_node, N)

    for (nidx, n) in enumerate(nodes)
        found = false
        for (zi, z) in enumerate(zones)
            if n ∈ zone_nodes[z]
                zone_of_node[nidx] = zi
                found = true
                break
            end
        end
        @assert found "Node $(n) not found in any zone"
    end


    zone_nodes_idx = [Int[] for _ in 1:Z]
    for nidx in 1:N
        z = zone_of_node[nidx]
        push!(zone_nodes_idx[z], nidx)
    end


    sets = Sets(H, Z, N, I, L, S, 1:H, 1:Z, 1:N, 1:I, 1:L, 1:S, zones, nodes, techs)
    maps = Maps(zone_of_node, zone_to_idx, tech_to_idx, zone_nodes_idx)
    return sets, maps
end


function define_parameters(data, load, pv, wind_on, nptdf, lines, weights, scarcity_df, sets::Sets, maps::Maps)
    H, Z, N, I, L = sets.H, sets.Z, sets.N, sets.I, sets.L

    # demand per zone
    D_zonal = Array{Float64}(undef, H, Z)
    for (zidx, z) in enumerate(sets.zones)
        D_zonal[:, zidx] = Float64.(safe_col(load, z))[1:H]
    end



    # elastic/price parameters per zone
    ela = zeros(Float64, Z)
    WTP = zeros(Float64, Z)
    for (zidx, z) in enumerate(sets.zones)
        ela[zidx] = Float64(get_any(data["Consumers"][z], "ela"))
        WTP[zidx] = Float64(get_any(data["Consumers"][z], "WTP"))
    end

    # demand share per zone-node (Z × N)
    S = zeros(Float64, Z, N)
    for (zidx, z) in enumerate(sets.zones)
        local_nodes = maps.zone_nodes_idx[zidx]
        if haskey(data["Consumers"][z], "NodeShare")
            raw = Dict(Symbol(k)=>Float64(v) for (k,v) in data["Consumers"][z]["NodeShare"])
            tot = sum(get(raw, sets.nodes[n], 0.0) for n in local_nodes)
            if tot == 0
                for n in local_nodes; S[zidx, n] = 1.0/length(local_nodes); end
            else
                for n in local_nodes; S[zidx, n] = get(raw, sets.nodes[n], 0.0) / tot; end
            end
        else
            for n in local_nodes; S[zidx, n] = 1.0/length(local_nodes); end
        end
        @assert isapprox(sum(S[zidx, local_nodes]), 1.0; atol=1e-10) "Demand shares for zone $(z) must sum to 1"
    end

    # availability factors H × I × N
    function zone_af(tech::Symbol, z::String)
        if tech in (:PV, :WindOnshore)
            df = tech === :PV ? pv : wind_on
            return Float64.(safe_col(df, z))[1:H]
        else
            return ones(Float64, H)
        end
    end
    
    AV = zeros(Float64, H, I, N)
    for (i, tech) in enumerate(sets.techs), (zidx, z) in enumerate(sets.zones)
        af = zone_af(tech, z)
        for nidx in maps.zone_nodes_idx[zidx]
            @inbounds AV[:, i, nidx] .= af
        end
    end


    # existing capacity per (i,n)
    y_node = zeros(Float64, I, N)
    for (zidx, z) in enumerate(sets.zones)
        local_nodes = maps.zone_nodes_idx[zidx]
        for (tname, genco) in data["Generators"][z]
            i = maps.tech_to_idx[Symbol(tname)]
            y0 = Float64(get_any(genco, "C"))
            if haskey(genco, "NodeShare")
                shares = Dict(Symbol(k)=>Float64(v) for (k,v) in genco["NodeShare"])
                total = sum(get(shares, sets.nodes[j], 0.0) for j in local_nodes)
                total = total == 0 ? 1.0 : total
                for n in local_nodes
                    s = get(shares, sets.nodes[n], 0.0) / total
                    y_node[i, n] += s * y0
                end
            else
                for n in local_nodes
                    y_node[i, n] += (1.0/length(local_nodes)) * y0
                end
            end
        end
    end


    # variable & investment costs per (i,z)
    MC = zeros(Float64, I, Z)
    A = zeros(Float64, I, Z)
    IC = zeros(Float64, I, Z)
    max_cap = zeros(Float64, I, Z)
    for (zidx, z) in enumerate(sets.zones), (i, tech) in enumerate(sets.techs)
        if haskey(data["Generators"][z], String(tech))
            genco = data["Generators"][z][String(tech)]
        else
            error("Tech $(tech) missing for zone $(z). Add it or provide a default; current code assumes every (i,z) exists.")
        end

        MC[i, zidx] = Float64(get_any(genco, "b"))
        A[i, zidx] = Float64(get_any(genco, "a"))
        IC[i, zidx] = Float64(get_any(genco, "I"))
        max_cap[i, zidx] = Float64(get_any(genco, "max_cap"))
    end

    # PTDF + limits
    PTDF = Matrix{Float64}(undef, L, N)
    for (j, n) in enumerate(sets.nodes)
        PTDF[:, j] = Float64.(safe_col(nptdf, String(n)))
    end
    # println(DataFrame(PTDF, :auto))

    Fmax = Float64.(safe_col(lines, :Fmax))
    @assert size(PTDF,1) == L "PTDF rows must match number of lines"

    margins = zeros(Float64, sets.Z)
    for (zidx, z) in enumerate(sets.zones)
        margins[zidx] = Float64(get_any(data["CM"][z], "capacity_margin"))
    end

    # Capacity demand per zone
    CD_ref = zeros(Float64, sets.Z)
    for (zidx, z) in enumerate(sets.zones)
        CD_ref[zidx] = Float64(get_any(data["CM"][z], "capacity_target"))
    end

    scarcity_matrix = Matrix{Float64}(scarcity_df[:, sets.zones])
    return Params(weights, D_zonal, AV, y_node, MC, max_cap, A, IC, S, PTDF, Fmax, margins, CD_ref, scarcity_matrix), ela, WTP
end

function build_planner_cm(; data, load, pv, wind_on, nodal_ptdf_df, lines, weights, scarcity_df)
    sets, maps = build_sets_and_maps(data, lines, scarcity_df)
    params, ela, WTP = define_parameters(data, load, pv, wind_on, nodal_ptdf_df, lines, weights, scarcity_df, sets, maps)

    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag"=>0))

    H, Z, N, I, L, S = sets.H, sets.Z, sets.N, sets.I, sets.L, sets.S
    JH, JZ, JN, JI, JL, JS = sets.JH, sets.JZ, sets.JN, sets.JI, sets.JL, sets.JS
    W, D_zonal, AV, y_node = params.W, params.D_zonal, params.AV, params.y_node
    MC, max_cap, A, IC, PTDF, Fmax, S = params.MC, params.max_cap, params.A, params.IC, params.PTDF, params.Fmax, params.S
    CD_ref, scarcity_matrix, margins = params.CD_ref, params.scarcity_matrix, params.margins
    zone_of_node = maps.zone_of_node
   
    # EOM Variables
    @variable(model, y[JI, JZ] >= 0)
    @variable(model, y_bar[JI, JN] >= 0)
    @variable(model, g[JH, JI, JZ] >= 0)
    @variable(model, g_bar[JH, JI, JN] >= 0)
    @variable(model, r[JH, JN])
    @variable(model, f[JH, JL])
    @variable(model, p[JH, JZ])
    @variable(model, d_inel[JH, JZ] >= 0)
    @variable(model, d_ela[JH, JZ]  >= 0)
    @variable(model, ens[JH, JZ] >= 0)

    ### CM variables
    @variable(model, cap_cm[JI, JZ] >= 0)           # capacity sold in CM per technology per zone
    @variable(model, cap_cm_bar[JI, JN] >= 0)       # capacity sold in CM at nodal level
    @variable(model, g_cm[JS, JI, JN] >= 0)         # capacity deployed at nodal level in scarcity
    @variable(model, r_cm[JS, JN])                  # nodal net position in CM scenarios
    @variable(model, f_cm[JS, JL])                  # line flows in CM scenarios
    @variable(model, p_cm[JZ])                      # zonal capacity net positions
    @variable(model, CD[JZ] >= 0)                   # capacity demand in CM per zone


    # negative utility function per zone
    function consumer_utility(t, z)
        w = WTP[z]
        e = ela[z]
        if e == 0.0 || D_zonal[t, z] == 0.0
            return -w * d_inel[t, z]
        else
            return -w * (d_inel[t, z] + d_ela[t, z]) + (w / (2 * e * D_zonal[t, z])) * d_ela[t, z]^2
        end
    end

 

    # marginal cost + quadratic cost
    @expression(model, gen_cost,
        sum(W[t,z] * MC[i, z] * g[t, i, z] for t in JH, i in JI, z in JZ) +
        sum(W[t,z] * A[i, z] / 2 * g[t, i, z]^2 for t in JH, i in JI, z in JZ)
    )
    # investment cost
    @expression(model, inv_cost, sum(IC[i, z] * y[i, z] for i in JI, z in JZ))

    # consumer utility
    @expression(model, neg_utility,  sum(W[t,z] * consumer_utility(t, z) for t in JH, z in JZ))

    # Objective: minimize (costs - utility)
    @objective(model, Min, (inv_cost + gen_cost + neg_utility))

    # Constraints

    # demand
    @constraint(model, [t in JH, z in JZ], d_inel[t, z] == (1 - ela[z]) * D_zonal[t, z] - ens[t, z])
    @constraint(model, [t in JH, z in JZ], d_ela[t, z]  <=  ela[z] * D_zonal[t, z])

    # zonal balance constraints
    bal = @constraint(model, bal[t in JH, z in JZ], -p[t, z] + sum(g[t, i, z] for i in JI) - (d_inel[t, z] + d_ela[t, z]) == 0) # imports negative

    # Aggregate nodal net positions to zonal net positions
    agg = @constraint(model, [t in JH, z in JZ], 
        p[t, z] == sum(r[t, n] for n in JN if zone_of_node[n] == z))

    # Nodal balance constraints
    nbal = @constraint(model, [t in JH, n in JN], 
        r[t, n] == sum(g_bar[t, i, n] for i in JI) - 
                   (d_inel[t, zone_of_node[n]] + d_ela[t, zone_of_node[n]]) * S[zone_of_node[n], n])

    # nodal capacity limits
    cap  = @constraint(model, [t in JH, i in JI, n in JN], g_bar[t, i, n] <= AV[t, i, n] * (y_node[i, n] + y_bar[i, n] ))

    # Mainly constrains renewables investment
    Ren_cap = @constraint(model, [i in JI, z in JZ], sum(y_node[i, n] for n in JN if zone_of_node[n] == z) + y[i, z] <= max_cap[i, z])

    # link between zonal and nodal capacity
    alloc = @constraint(model, [i in JI, z in JZ], y[i, z] == sum(y_bar[i, n] for n in JN if zone_of_node[n] == z))

    # DC power flow constraints
    fmap  = @constraint(model, [t in JH, l in JL], f[t, l] == sum(PTDF[l, n] * r[t, n] for n in JN))

    # thermal limits
    therm = @constraint(model, [t in JH, l in JL], -Fmax[l] <= f[t,l] <= Fmax[l])

    # system balance
    sbal  = @constraint(model, [t in JH], sum(r[t, n] for n in JN) == 0)

    # zonal capacity limits
    zonal_cap = @constraint(model, [t in JH, i in JI, z in JZ], g[t, i, z] <= AV[t, i, z] * (y[i, z] + sum(y_node[i, n] for n in JN if zone_of_node[n] == z)))



    # initialize
    cm_bal = nothing
    cm_cap = nothing
    cm_alloc = nothing
    cm_gen = nothing
    cm_nbal = nothing
    cm_agg = nothing
    cm_fmap = nothing
    cm_therm = nothing
    cm_sbal = nothing
    cm_gbal = nothing
    cm_netpos = nothing
    cm_atc_limit = nothing


    ### capacity market constraints (Flow-based)

    if data["Network"]["coupling"] == "FB"

        @variable(model, sys_req >= 0) # newly added

        # Zonal balance in CM: capacity offered + imports = capacity demand
        cm_bal = @constraint(model, [z in JZ],  sum(cap_cm[i, z] for i in JI) + p_cm[z] == CD[z])

        # # Capacity limits:
        cm_cap = @constraint(model, [i in JI, n in JN], cap_cm_bar[i, n] <= y_bar[i, n])
        
        # Link between zonal and nodal capacity in CM -> sum of nodal capacity offer in CM == nodal capacity sold in CM
        cm_alloc = @constraint(model, [i in JI, z in JZ], cap_cm[i, z] == sum(cap_cm_bar[i, n] for n in JN if zone_of_node[n] == z))

        # nodal capacity generation in CM <= derated capacity offer in CM
        # cm_gen = @constraint(model, [s in JS, i in JI, n in JN], g_cm[s, i, n] <= cap_cm_bar[i, n])
        cm_gen = @constraint(model, [s in JS, i in JI, n in JN], g_cm[s, i, n] <= y_bar[i, n])


        nodal_capacity_demand = @expression(model, [s in JS, n in JN], S[zone_of_node[n], n] * CD[zone_of_node[n]] * scarcity_matrix[s, zone_of_node[n]])
        nodal_capacity_demand_ref = @expression(model, [s in JS, n in JN], S[zone_of_node[n], n] * CD_ref[zone_of_node[n]] * scarcity_matrix[s, zone_of_node[n]])
            
        zonal_req = @expression(model, [s in JS, z in JZ], sum(nodal_capacity_demand_ref[s, n] for n in JN if zone_of_node[n] == z))
        
        # System requirement equals total CD across zones
        cm_demand = @constraint(model, sum(CD[z] for z in JZ) == sys_req)
        
        # System requirement must be sufficient for each scenario
        cm_adequacy = @constraint(model, [s in JS], sys_req >= sum(zonal_req[s, z] for z in JZ))

        # nodal balance
        cm_nbal = @constraint(model, [s in JS, n in JN], r_cm[s, n] == sum(g_cm[s, i, n] for i in JI) - nodal_capacity_demand[s, n])

        # capacity net position in CM = sum of nodal capacity net positions in CM
        cm_agg = @constraint(model, [s in JS, z in JZ], p_cm[z] >= - sum(r_cm[s, n] for n in JN if zone_of_node[n] == z))

        # DC power flow for power delivery during scarcity scenario in CM
        cm_fmap = @constraint(model, [s in JS, l in JL], f_cm[s, l] == sum(PTDF[l, n] * r_cm[s, n] for n in JN))

        # Thermal limits in scarcity scenarios
        cm_therm = @constraint(model, [s in JS, l in JL], -Fmax[l] <= f_cm[s, l] <= Fmax[l])

        # System balance in scarcity scenarios -> sum of nodal net positions in CM == 0
        cm_sbal = @constraint(model, [s in JS], sum(r_cm[s, n] for n in JN) == 0)
        
        cm_gbal = @constraint(model, sum(p_cm[z] for z in JZ) == 0)

    # Capacity market constraints (ATCMC)
    elseif data["Network"]["coupling"] == "ATC"
        # Parse ATC data from config
        atc_data = data["Network"]["ATC"]
        ATC = Dict{Tuple{Symbol,Symbol}, Tuple{Float64,Float64}}()
        
        # Parse TCONNECT data - interconnected zones
        TCONNECT = [(Symbol(t[1]), Symbol(t[2])) for t in data["Network"]["TCONNECT"]]
        
        # Build ATC dictionary with forward/backward limits
        for from_zone in keys(atc_data)
            for (to_zone, limits) in atc_data[from_zone]
                from_sym = Symbol(from_zone)
                to_sym = Symbol(to_zone)
                forward_limit = Float64(limits[1])
                backward_limit = Float64(limits[2])
                ATC[(from_sym, to_sym)] = (forward_limit, backward_limit)
            end
        end
        
        # Create capacity exchange variables
        @variable(model, ex_cm[t in TCONNECT], base_name="ex_cm")
        
        # Zone to symbol mapping for constraint reference
        zone_syms = Symbol.(sets.zones)
                
        # zonal balance in CM -> capacity sold in CM per zone = capacity manager net position in CM + capacity demand in zone
        cm_bal = @constraint(model, [z in JZ], sum(cap_cm[i, z] for i in JI) - CD[z] + p_cm[z] == 0) # imports positive

        # Capacity limits:
        cm_cap = @constraint(model, [i in JI, n in JN], cap_cm_bar[i, n] <= y_bar[i, n])
        
        # Link between zonal and nodal capacity in CM -> sum of nodal capacity offer in CM == nodal capacity sold in CM
        cm_alloc = @constraint(model, [i in JI, z in JZ], cap_cm[i, z] == sum(cap_cm_bar[i, n] for n in JN if zone_of_node[n] == z))

        # nodal capacity generation in CM <= derated capacity offer in CM
        # cm_gen = @constraint(model, [s in JS, i in JI, n in JN], g_cm[s, i, n] <= cap_cm_bar[i, n])
        cm_gen = @constraint(model, [s in JS, i in JI, n in JN], g_cm[s, i, n] <= y_bar[i, n])

        # Global balance constraint (sum of net positions = 0)
        cm_gbal = @constraint(model, sum(p_cm[z] for z in JZ) == 0)
        
        # Define zonal net positions based on exchanges
        cm_netpos = @constraint(model, [z in JZ],
            p_cm[z] == 
            sum(ex_cm[t] for t in TCONNECT if t[2] == zone_syms[z]) -
            sum(ex_cm[t] for t in TCONNECT if t[1] == zone_syms[z])
        )
        
        # ATC limits on exchanges
        cm_atc_limit = @constraint(model, [t in TCONNECT], 
            ATC[t][2] <= ex_cm[t] <= ATC[t][1]
        )

        cm_demand_lower = @constraint(model, [z in JZ], CD[z] >= CD_ref[z])
        cm_demand_upper = @constraint(model, [z in JZ], CD[z] <= CD_ref[z] * (1 + margins[z]))
        
    else
        # constrain capacity trade to zero
        @constraint(model, [z in JZ], p_cm[z] == 0)
                # zonal balance in CM -> capacity sold in CM per zone = capacity manager net position in CM + capacity demand in zone
        cm_bal = @constraint(model, [z in JZ], sum(cap_cm[i, z] for i in JI) - CD[z] == 0) # imports positive

        # cap_cm must be less than installed capacity in the zone
        cm_cap = @constraint(model, [i in JI, z in JZ], cap_cm[i, z] <= sum(y_bar[i, n] for n in JN if zone_of_node[n] == z))

        cm_demand_lower = @constraint(model, [z in JZ], CD[z] >= CD_ref[z])
        cm_demand_upper = @constraint(model, [z in JZ], CD[z] <= CD_ref[z] * (1 + margins[z]))
    end



    model.ext = Dict{Symbol,Any}(
        :sets => sets,
        :maps => maps,
        :params => params,
        :vars => Dict(
            :y => y, :y_bar => y_bar, :g => g, :g_bar => g_bar,
            :r => r, :f => f, :p => p, :d_inel => d_inel, :d_ela => d_ela,
            :cap_cm => cap_cm, :cap_cm_bar => cap_cm_bar, :g_cm => g_cm, 
            :r_cm => r_cm, :f_cm => f_cm, :p_cm => p_cm, :CD => CD
        ),
        :constraint => Dict(
            :bal => bal, :agg => agg, :nbal => nbal, :cap => cap,
            :alloc => alloc, :fmap => fmap, :therm => therm, :sbal => sbal,
            :zonal_cap => zonal_cap, :Ren_cap => Ren_cap,
            :cm_bal => cm_bal, 
            :cm_cap => cm_cap, :cm_alloc => cm_alloc,
            :cm_gen => cm_gen, :cm_nbal => cm_nbal, :cm_agg => cm_agg,
            :cm_fmap => cm_fmap, :cm_therm => cm_therm, 
            :cm_sbal => cm_sbal, :cm_gbal => cm_gbal,
        ),
        :scalars => Dict(:ela => ela, :WTP => WTP),
        :expressions => Dict(:ens => ens, :gen_cost => gen_cost, 
            :inv_cost => inv_cost, :neg_utility => neg_utility)
    )
    return model
end

function solve_and_save(
    data::Dict,
    load::DataFrame,
    pv::DataFrame,
    wind_on::DataFrame,
    nptdf::DataFrame,
    lines::DataFrame,
    weights,
    scarcity_df::DataFrame;
    output_dir = joinpath(@__DIR__, "Results", "Planner_CM"),
)
    mkpath(output_dir)

    # Build and solve
    model = build_planner_cm(;
        data=data,
        load=load,
        pv=pv,
        wind_on=wind_on,
        nodal_ptdf_df=nptdf,
        lines=lines,
        weights=weights,
        scarcity_df=scarcity_df,
    )
    optimize!(model)

    term = termination_status(model)
    if term != MOI.OPTIMAL
        ps = primal_status(model)
        error("Solver did not return OPTIMAL. termination_status=$(term), primal_status=$(ps)")
    end

    ext    = model.ext
    sets   = ext[:sets]
    maps   = ext[:maps]
    params = ext[:params]
    expressions = ext[:expressions]

    JH, JI, JZ, JN, JL, JS = sets.JH, sets.JI, sets.JZ, sets.JN, sets.JL, sets.JS
    zones, nodes, techs = sets.zones, sets.nodes, sets.techs

    vars = ext[:vars]
    constraints = ext[:constraint]
    @assert all(haskey(vars, s) for s in (:y, :y_bar, :g_bar, :p, :f, :d_inel, :d_ela, :cap_cm, :cap_cm_bar, :g_cm, :p_cm)) "model.ext.vars is missing expected fields"

    y_var     = vars[:y]        # (I,Z)   installed zonal capacity
    ybar_var  = vars[:y_bar]    # (I,N)   nodal allocation of new capacity
    g_var     = vars[:g]        # (H,I,N) nodal generation
    gbar_var  = vars[:g_bar]    # (H,I,N) nodal dispatch
    p_var     = vars[:p]        # (H,Z)   zonal net positions
    f_var     = vars[:f]        # (H,L)   line flows
    dinel_var = vars[:d_inel]   # (H,Z)   inelastic demand served
    dela_var  = vars[:d_ela]    # (H,Z)   elastic demand served
    capcm_var = vars[:cap_cm]      # (I,Z)   capacity sold in CM per technology per zone
    capcmbar_var = vars[:cap_cm_bar]  # (I,N)   capacity sold in CM at nodal level
    gcm_var   = vars[:g_cm]        # (S,I,N) capacity deployed at nodal level in scarcity
    pcm_var   = vars[:p_cm]        # (Z)     zonal capacity net positions
    CD_var   = vars[:CD]        # (Z)     capacity demand in CM per zone


    capacity     = value.(y_var)
    aux_capacity  = value.(ybar_var)
    gen          = value.(g_var)
    aux_dispatch = value.(gbar_var)
    net_pos = value.(p_var)
    flow     = value.(f_var)
    inel   = value.(dinel_var)
    elastic  = value.(dela_var)
    cap_cm   = value.(capcm_var)
    cap_cm_bar = value.(capcmbar_var)
    # println(cap_cm_bar)
    gen_cm   = value.(gcm_var)
    pos_cm   = value.(pcm_var)
    capacity_demand = value.(CD_var)

    # Zonal prices (€/MWh) from balance duals, unweighted
    W = params.W
    rho = [dual(constraints[:bal][t, z]) / W[t,z] for t in JH, z in JZ]
    cm_price = [dual(constraints[:cm_bal][z]) for z in JZ]
    # extract dual of alloc
    dual_alloc = [dual(constraints[:alloc][i, z]) for i in JI, z in JZ] # show the dataframe in the terminal
    println(DataFrame(dual_alloc, :auto))

    # Existing capacity C[i,z] from input data (for totals)
    Cz = zeros(length(JI), length(JZ))
    for (zname, zidx) in maps.zone_to_idx
        for (tname, pars) in data["Generators"][zname]
            i = maps.tech_to_idx[Symbol(tname)]
            Cz[i, zidx] = Float64(get_any(pars, "C"))
        end
    end

    # Calculate total capacity offered per zone
    total_cap_offer = [sum(cap_cm[i, z] for i in JI) for z in JZ]
    capacity_manager = [pos_cm[z] for z in JZ] # Zonal capacity net positions
    

    for z in zones
        zidx = maps.zone_to_idx[z]
        df = DataFrame(Timestep = collect(JH))
        df.EOM_price = rho[:, zidx]
        df[!, :CM_price] = fill(cm_price[zidx], length(JH))
        df[!, :TotalCapOffer] = fill(total_cap_offer[zidx], length(JH))

            # Get zonal generation per tech
            for tech in techs
                i = maps.tech_to_idx[tech]
                gen_t = gen[:, i, zidx]  # Direct access to zonal generation
                tstr = String(tech)
                df[!, Symbol("Gen_$(z)_$(tstr)")]                = gen_t
                df[!, Symbol("new_Capacity_Gen_$(z)_$(tstr)")]   = fill(capacity[i, zidx], length(JH))
                df[!, Symbol("Capacity_Gen_$(z)_$(tstr)")] = fill(Cz[i, zidx] + capacity[i, zidx], length(JH))
            end

        cons_total = inel[:, zidx] .+ elastic[:, zidx]

        ens = value.(expressions[:ens])[:, zidx]

        df[!, Symbol("Cons_$(z)")]           = cons_total
        df[!, Symbol("Inelastic_Cons_$(z)")] = inel[:, zidx]
        df[!, Symbol("Elastic_Cons_$(z)")]   = elastic[:, zidx]
        df[!, Symbol("CapDemand_Cons_$(z)")] = fill(capacity_demand[zidx], length(JH))
        df[!, :NetworkManager]               = -1 * net_pos[:, zidx] # imports positive
        df[!, :CapacityManager] = fill(capacity_manager[zidx], length(JH))
        df[!, Symbol("ENS_Cons_$(z)")]       = ens


        for col in names(df)
            if eltype(df[!, col]) <: Number
                df[!, col] = round.(df[!, col], digits=2)
            end
        end

        CSV.write(joinpath(output_dir, "planner_zone_$(z).csv"), df, delim=";")
    end

    price_df = DataFrame(Timestep = collect(JH))
    for (zidx, z) in enumerate(zones)
        price_df[!, Symbol("Price_$(z)")] = rho[:, zidx]
    end

    invest_df = DataFrame(Technology = String.(techs))
    for (zidx, z) in enumerate(zones)
        invest_df[!, Symbol("Invest_$(z)")] = capacity[:, zidx]
    end

    nodal_invest_df = DataFrame(
        Technology = repeat(String.(techs), sets.N),
        Node       = repeat(String.(String.(nodes)), inner = sets.I),
        Allocation = vec(aux_capacity),
    )

    dispatch_df = DataFrame(
        Timestep   = repeat(collect(JH), outer = sets.I * sets.N),
        Technology = repeat(repeat(String.(techs), inner = sets.H), outer = sets.N),
        Node       = repeat(String.(String.(nodes)), inner = sets.H * sets.I),
        Generation = vec(aux_dispatch),
    )

    cm_results_df = DataFrame(
        Zone = zones,
        CM_price = [cm_price[maps.zone_to_idx[z]] for z in zones],
        TotalCapOffer = [total_cap_offer[maps.zone_to_idx[z]] for z in zones],
        CapacityManager = [capacity_manager[maps.zone_to_idx[z]] for z in zones]
    )


    flow_df = DataFrame(Timestep = repeat(collect(JH), outer = sets.L))
    line_ids = :line_id ∈ names(lines) ? lines[!, :line_id] : 1:nrow(lines)
    flow_df[!, :Line] = repeat(line_ids, inner = sets.H)
    flow_df[!, :Flow] = vec(flow)

    return Dict(
        :status               => term,
        :objective            => objective_value(model),
        :prices_df            => price_df,
        :investments_df       => invest_df,
        :nodal_investments_df => nodal_invest_df,
        :dispatch_df          => dispatch_df,
        :flows_df             => flow_df,
        :cm_results_df        => cm_results_df,
        :out_dir              => output_dir,
        :model                => model,
    )
end

inputs = read_inputs(@__DIR__)
cp_cm = solve_and_save(
    inputs.data, 
    inputs.load, 
    inputs.pv, 
    inputs.wind_on, 
    inputs.nptdf, 
    inputs.lines, 
    inputs.weights,
    inputs.scarcity;
    output_dir=joinpath(@__DIR__, "Results", "Planner_CM")
)
@show objective_value(cp_cm[:model])

# Get y_node (existing capacity parameter)
model = cp_cm[:model]
params = model.ext[:params]
sets = model.ext[:sets]
maps = model.ext[:maps]
y_node = params.y_node  # This is the existing capacity matrix (I×N)

# Get y_bar (new capacity allocation variable)
vars = model.ext[:vars]
y_bar_var = vars[:y_bar]
y_bar_values = value.(y_bar_var)  # This is the optimal new capacity allocation (I×N)

# Create DataFrame for existing capacity
sets = cp_cm[:model].ext[:sets]
techs = sets.techs
nodes = sets.nodes

existing_capacity_df = DataFrame(
    Technology = repeat(String.(techs), sets.N),
    Node = repeat(String.(nodes), inner = sets.I),
    ExistingCapacity = vec(y_node)
)

# Create DataFrame for new capacity allocation
new_capacity_df = DataFrame(
    Technology = repeat(String.(techs), sets.N),
    Node = repeat(String.(nodes), inner = sets.I),
    NewCapacity = vec(y_bar_values)
)

# Add y_node to y_bar to get total capacity per (i,n)
total_capacity_df = copy(existing_capacity_df)
total_capacity_df[!, :TotalCapacity] = existing_capacity_df[!, :ExistingCapacity] .+ new_capacity_df[!, :NewCapacity]

# Calculate total capacity per node (summing across all technologies)
node_summary_df = combine(
    groupby(total_capacity_df, :Node),
    :ExistingCapacity => sum => :ExistingCapacity,
    :TotalCapacity => sum => :TotalCapacity
)

# Calculate new capacity per node
node_summary_df[!, :NewCapacity] = node_summary_df[!, :TotalCapacity] .- node_summary_df[!, :ExistingCapacity]

# Print results
# println("Total capacity by technology and node:")
# println(total_capacity_df)
println("\nTotal capacity by node:")
println(node_summary_df)
