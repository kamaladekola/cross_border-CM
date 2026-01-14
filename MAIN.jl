## 0. Set-up code
# Home directory
const home_dir = @__DIR__

######### Add packages #########
# import Pkg
# Pkg.add.(["DataStructures","ProgressBars","TimerOutputs","ArgParse","JuMP","Gurobi","CSV","YAML","DataFrames","JLD2"])


# Include packages 
using JuMP, Gurobi
using DataFrames, CSV, YAML, DataStructures
using ProgressBars, Printf
using TimerOutputs # profiling 
using Base.Threads: @spawn 
using Base: split
using ArgParse # Parsing arguments from the command line

# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")   
GRBsetparam(GUROBI_ENV, "Threads", "4")
println("        ")

# Include functions
include(joinpath(home_dir,"Source","utility_functions.jl"))
include(joinpath(home_dir,"Source","define_common_parameters.jl"))
include(joinpath(home_dir,"Source","define_EOM_parameters.jl"))
include(joinpath(home_dir,"Source","define_consumer_parameters.jl"))
include(joinpath(home_dir,"Source","define_generator_parameters.jl"))
include(joinpath(home_dir,"Source","define_interconnector_parameters.jl")) # added interconnector parameters
include(joinpath(home_dir,"Source","define_capacityIC_parameters.jl")) # added capacity manager parameters
include(joinpath(home_dir,"Source","build_consumer_agent.jl"))
include(joinpath(home_dir,"Source","build_generator_agent.jl"))
include(joinpath(home_dir,"Source","build_interconnector_agent.jl")) # added interconnector agent
include(joinpath(home_dir,"Source","build_capacityIC_agent.jl")) # added capacity manager agent
include(joinpath(home_dir,"Source","define_results.jl"))
include(joinpath(home_dir,"Source","ADMM.jl"))
include(joinpath(home_dir,"Source","ADMM_subroutine.jl"))
include(joinpath(home_dir,"Source","solve_consumer_agent.jl"))
include(joinpath(home_dir,"Source","solve_generator_agent.jl"))
include(joinpath(home_dir,"Source","solve_interconnector_agent.jl")) # added interconnector agent
include(joinpath(home_dir,"Source","solve_capacityIC_agent.jl")) # added capacity manager agent
include(joinpath(home_dir,"Source","define_getATC.jl"))
include(joinpath(home_dir,"Source","build_getATC.jl"))
include(joinpath(home_dir,"Source","solve_getATC.jl"))
include(joinpath(home_dir,"Source","update_rho.jl"))
include(joinpath(home_dir,"Source","save_results.jl"))

# Data common to all scenarios data 
data = YAML.load_file(joinpath(home_dir,"Input","config.yaml"))
weights = CSV.read(joinpath(home_dir,"Input","weights.csv"),delim=";",DataFrame) # to be removed
load = CSV.read(joinpath(home_dir,"Input","load.csv"),delim=";",DataFrame) # columns are zones
pv = CSV.read(joinpath(home_dir,"Input","pv.csv"),delim=";",DataFrame)
# wind_offshore = CSV.read(joinpath(home_dir,"Input","wind_offshore.csv"),delim=";",DataFrame)
wind_onshore = CSV.read(joinpath(home_dir,"Input","wind_onshore.csv"),delim=";",DataFrame)
ptdf = CSV.read(joinpath(home_dir,"Input","ptdf.csv"),delim=";",DataFrame)
nodal_ptdf = CSV.read(joinpath(home_dir,"Input","nodal_ptdf.csv"),delim=";",DataFrame)
lines = CSV.read(joinpath(home_dir,"Input","lines.csv"),delim=";",DataFrame)
participation_matrix = CSV.read(joinpath(home_dir,"Input","participation_matrix.csv"),delim=";",DataFrame)
derating_factor = CSV.read(joinpath(home_dir,"Input","derating_factor.csv"),delim=";",DataFrame)
scarcity = CSV.read(joinpath(home_dir,"Input","scarcity.csv"),delim=";",DataFrame)

# Overview scenarios
scenario_overview = CSV.read(joinpath(home_dir,"overview_scenarios.csv"),DataFrame,delim=";")
sensitivity_overview = CSV.read(joinpath(home_dir,"overview_sensitivity.csv"),DataFrame,delim=";") 

# Create file with results 
# add column for sensitivity analsysis
if isfile(joinpath(home_dir,string("overview_results.csv"))) != 1
    CSV.write(joinpath(home_dir,string("overview_results.csv")),DataFrame(),delim=";",header=["scen_number";"sensitivity";"n_iter";"walltime";"PrimalResidual_EOM"; "DualResidual_EOM"])
end

# Create folder for results
if isdir(joinpath(home_dir,string("Results"))) != 1
    mkdir(joinpath(home_dir,string("Results")))
end

scen_number = 1 # for debugging purposes, comment the for-loop and replace it by a explicit definition of the scenario you'd like to study
# for scen_number in range(start_scen,stop=stop_scen,step=1)

println("    ")
println(string("######################                  Scenario ",scen_number,"                 #########################"))

## 1. Read associated input for this simulation
scenario_overview_row = scenario_overview[scen_number,:]
data = YAML.load_file(joinpath(home_dir,"Input","config.yaml")) # reload data to avoid previous sensitivity analysis affected data

if scenario_overview_row["Sens_analysis"] == "YES"  
    numb_of_sens = length((sensitivity_overview[!,:Parameter]))
else
    numb_of_sens = 0 
end    
sens_number = 1 # for debugging purposes, comment the for-loop and replace it by a explicit definition of the sensitivity you'd like to study
# for sens_number in range(1,stop=numb_of_sens+1,step=1) 
if sens_number >= 2
    println("    ") 
    println(string("#                                  Sensitivity ",sens_number-1,"                                      #"))
    parameter = split(sensitivity_overview[sens_number-1,:Parameter])
    if length(parameter) == 2
        data[parameter[1]][parameter[2]] = sensitivity_overview[sens_number-1,:Scaling]*data[parameter[1]][parameter[2]]
    elseif length(parameter) == 3
        data[parameter[1]][parameter[2]][parameter[3]] = sensitivity_overview[sens_number-1,:Scaling]*data[parameter[1]][parameter[2]][parameter[3]]
    else
        printnl("warning! Sensitivity analysis is not well defined!")
    end
end

println("    ")
println("Including all required input data: done")
println("   ")

## 2. Initiate models for representative agents

# Create an ordered list of zones to ensure consistent ordering throughout the code
zones = sort(string.(collect(keys(data["Consumers"])))) # zones = ["A","B","C"]

# Parameters/variables EOM
EOM = Dict()
CM = Dict()

EOM["nZones"] = length(zones)

println("Existing zones:" , zones)
println("   ")

agents = Dict{Symbol, Any}()

agents[:Gen] = String[]
agents[:Cons] = String[]
agents[:IC] = String[]

agents[:Gen] = ["Gen_$(z)_$(tech)" for z in zones for tech in keys(data["Generators"][z])] # generator agents for each zone and technology
agents[:Cons] = ["Cons_$(z)" for z in zones] # consumer agents for each zone
agents[:IC] = ["NetworkManager"] # interconnectors
agents[:CIC] = ["CapacityManager"] # Available capacity under scarcity


agents[:all] = union(agents[:Gen],agents[:Cons], agents[:IC], agents[:CIC])                                         # all agents in the game  
agents[:eom] = union(agents[:Gen],agents[:Cons], agents[:IC])                                         # agents participating in the EOM                           
agents[:cm] = union(agents[:Gen], agents[:Cons], agents[:CIC])                                                                     # agents participating in the CM -> may add agents[:Cons] if demand response participates in the CM

# grouping agents by zone
agents[:Gen_Z] = Dict(z => [m for m in agents[:Gen] if parse_agent_name(m)[1] == z] for z in zones)
agents[:Cons_Z] = Dict(z => [m for m in agents[:Cons] if parse_agent_name(m)[1] == z] for z in zones)


agents[:zones] = Dict(z => union(agents[:Gen_Z][z], agents[:Cons_Z][z]) for z in zones)       # agents in zone Z
agents[:eom_Z] = Dict(z => union(agents[:Gen_Z][z], agents[:Cons_Z][z]) for z in zones)       # agents participating in zone Z participating in the EOM
agents[:cm_Z] = Dict(z => union(agents[:Gen_Z][z], agents[:Cons_Z][z]) for z in zones)                                              # agents participating in zone Z participating in the CM -> may add agents[:Cons_Z] if demand response participates in the CM

# create one model per agent
mdict = Dict{String, Model}(i => Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))) for i in agents[:all])

# consumer models
for m in agents[:Cons]
    zone, _ = parse_agent_name(m)
    cons_data = merge(data["General"], data["Consumers"][zone], data["CM"][zone])

    define_common_parameters!(m, mdict[m], data, agents, scenario_overview_row, zones, ptdf, nodal_ptdf, lines, participation_matrix, derating_factor) # Parameters common to all agents
    define_consumer_parameters!(mdict[m], cons_data, load, weights)                            # Consumers
end

# Generator models
for m in agents[:Gen]
    zone, tech = parse_agent_name(m)
    if tech == "WindOnshore"
        af = wind_onshore
    elseif tech == "PV"
        af = pv
    else
        af = DataFrame(A = ones(data["General"]["nTimesteps"]), 
        B = ones(data["General"]["nTimesteps"]), 
        C = ones(data["General"]["nTimesteps"]))
    end
    gen_data = merge(data["General"], data["Generators"][zone][tech])
    define_common_parameters!(m, mdict[m], data, agents, scenario_overview_row, zones, ptdf, nodal_ptdf, lines, participation_matrix, derating_factor) # Parameters common to all agents
    define_generator_parameters!(mdict[m], gen_data, weights, af, zone)                             # Generators
end

# Interconnector models
for m in agents[:IC]
    IC_data = merge(data["General"], data["Network"], data["Consumers"])
    define_common_parameters!(m, mdict[m], data, agents, scenario_overview_row, zones, ptdf, nodal_ptdf, lines, participation_matrix, derating_factor) # Parameters common to all agents
    define_interconnector_parameters!(mdict[m], IC_data, zones, ptdf, nodal_ptdf, lines, weights)                # Interconnectors
end

# Capacity manager models
for m in agents[:CIC]
    define_common_parameters!(m, mdict[m], data, agents, scenario_overview_row, zones, ptdf, nodal_ptdf, lines, participation_matrix, derating_factor) # Parameters common to all agents
    define_capacityIC_parameters!(mdict[m], data, zones, ptdf, scarcity) # Capacity manager
    if data["Network"]["coupling"] == "ATC"
            define_getATC!(mdict[m]) # ATC parameters
    end
end

## 3. Define parameters for markets and representative agents

# define_EOM_parameters!(EOM,data,load,scenario_overview_row,zones)

# Calculate number of agents in each market
EOM["nAgents"] = length(agents[:eom])
EOM["nAgents_z"] = Dict(z => length(agents[:eom_Z][z]) for z in zones)
CM["nAgents"] = length(agents[:cm])
CM["nAgents_z"] = Dict(z => length(agents[:cm_Z][z]) for z in zones)
# println("Number of agents per zone: ", EOM["nAgents_z"])


println("Inititate model, sets and parameters: done")
println("   ")

## 4. Build models
for m in agents[:Cons]
    build_consumer_agent!(mdict[m], m, zones)
end
for m in agents[:Gen]
    build_generator_agent!(mdict[m], m, zones)
end
for m in agents[:IC]
    build_interconnector_agent!(mdict[m])
    
end
for m in agents[:CIC]
    build_capacityIC_agent!(mdict[m])
    # if data["Network"]["coupling"] == "ATC"
    #     build_getATC!(mdict[m])
    # end   
end

println("Build model: done")
println("   ")

## 5. ADMM process to calculate equilibrium
println("Find equilibrium solution...")
println("   ")
println("(Progress indicators on primal residuals, relative to tolerance: <1 indicates convergence)")
println("   ")

results = Dict()
ADMM = Dict()
TO = TimerOutput()
define_results!(merge(data["General"],data["ADMM"]),results,ADMM,agents,zones)           # initialize structure of results, only those that will be stored in each iteration
ADMM!(results,ADMM,EOM,CM,mdict,agents,scenario_overview_row,data,TO,zones)                 # calculate equilibrium 
ADMM["walltime"] =  TimerOutputs.tottime(TO)*10^-9/60                              # wall time 

println(string("Done!"))
println(string("        "))
println(string("Required iterations: ",ADMM["n_iter"]))
println(string("        "))

# per zone?
# println(string("RP EOM: ",  ADMM["Residuals"]["Primal"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
# println(string("RD EOM: ",  ADMM["Residuals"]["Dual"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
for zone in zones
    println(string("RP EOM for zone ", zone, ": ", ADMM["Residuals"]["Primal"]["EOM"][zone][end], " -- Tolerance: ", ADMM["Tolerance"]["EOM"]))
    println(string("RD EOM for zone ", zone, ": ", ADMM["Residuals"]["Dual"]["EOM"][zone][end], " -- Tolerance: ", ADMM["Tolerance"]["EOM"]))
    println(string("RP CM for zone ", zone, ": ", ADMM["Residuals"]["Primal"]["CM"][zone][end], " -- Tolerance: ", ADMM["Tolerance"]["CM"]))
    println(string("RD CM for zone ", zone, ": ", ADMM["Residuals"]["Dual"]["CM"][zone][end], " -- Tolerance: ", ADMM["Tolerance"]["CM"]))
end

println(string("        "))

## 6. Postprocessing and save results 
if sens_number >= 2
save_results(mdict,EOM,ADMM,results,data,agents,scenario_overview_row,sensitivity_overview[sens_number-1,:remarks], zones) 
# @save joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_",sensitivity_overview[sens_number-1,:remarks]))
else
save_results(mdict,EOM,ADMM,results,data,agents,scenario_overview_row,"ref", zones) 
# @save joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_ref"))
end

println("Postprocessing & save results: done")
println("   ")

# end # end loop over sensititivity
# end # end for loop over scenarios

println(string("##############################################################################################"))

