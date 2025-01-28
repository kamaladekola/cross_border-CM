## Topic: Go-E
# Author: Kenneth Bruninx
# Last update: November 2022

## 0. Set-up code
# HPC or not?
HPC = "NA" # NA (not applicable) or DelftBlue  

# Home directory
const home_dir = @__DIR__

if HPC == "DelftBlue"  # only for running this on DelftBlue
    ENV["GRB_LICENSE_FILE"] = "./Hpc/gurobi.lic"
    ENV["GUROBI_HOME"] = "./scratch/kbruninx/gurobi950/linux64"
    println(string("Number of threads: ", Threads.nthreads()))
end

if HPC == "ThinKing"  # only for running this on VSC
    # ENV["GRB_LICENSE_FILE"] = " "
    # ENV["GUROBI_HOME"] = " "
end

######### Add packages #########
import Pkg
# Pkg.add("DataStructures")
# Pkg.add("ProgressBars")
# Pkg.add("TimerOutputs")
Pkg.add("ArgParse")


# Include packages 
using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using TimerOutputs # profiling 
using Base.Threads: @spawn 
using Base: split
using ArgParse # Parsing arguments from the command line
# using JLD2 # save workspace

# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")   
GRBsetparam(GUROBI_ENV, "Threads", "4")   
println("        ")

# Include functions
include(joinpath(home_dir,"Source","define_common_parameters.jl"))
include(joinpath(home_dir,"Source","define_EOM_parameters.jl"))
include(joinpath(home_dir,"Source","define_consumer_parameters.jl"))
include(joinpath(home_dir,"Source","define_generator_parameters.jl"))
include(joinpath(home_dir,"Source","build_consumer_agent.jl"))
include(joinpath(home_dir,"Source","build_generator_agent.jl"))
include(joinpath(home_dir,"Source","define_results.jl"))
include(joinpath(home_dir,"Source","ADMM.jl"))
include(joinpath(home_dir,"Source","ADMM_subroutine.jl"))
include(joinpath(home_dir,"Source","solve_consumer_agent.jl"))
include(joinpath(home_dir,"Source","solve_generator_agent.jl"))
include(joinpath(home_dir,"Source","update_rho.jl"))
include(joinpath(home_dir,"Source","save_results.jl"))

# Data common to all scenarios data 
data = YAML.load_file(joinpath(home_dir,"Input","config.yaml"))
ts = CSV.read(joinpath(home_dir,"Input","timeseries.csv"),delim=";",DataFrame)

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

# Scenario number 
if HPC == "DelftBlue"  
   function parse_commandline()
       s = ArgParseSettings()
       @add_arg_table! s begin
           "--start_scen"
               help = "Enter the number of the first scenario here"
               arg_type = Int
               default = 1
            "--stop_scen"
               help = "Enter the number of the last scenario here"
               arg_type = Int
               default = 1
       end
       return parse_args(s)
   end
   # Simulation number as argument:
   dict_sim_number =  parse_commandline()
   start_scen = dict_sim_number["start_scen"]
   stop_scen = dict_sim_number["stop_scen"]
else
    # Range of scenarios to be simulated
    start_scen = 1
    stop_scen = 2
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
agents = Dict()
agents[:Gen] = [id for id in keys(data["Generators"])] 
agents[:Cons] = [id for id in keys(data["Consumers"])]
agents[:all] = union(agents[:Gen],agents[:Cons]) # all agents in the game  
agents[:eom] = union(agents[:Gen],agents[:Cons]) # agents participating in the EOM                           
mdict = Dict(i => Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))) for i in agents[:all])

## 3. Define parameters for markets and representative agents
# Parameters/variables EOM
EOM = Dict()
define_EOM_parameters!(EOM,data,ts,scenario_overview_row)
# Market 2

# consumer models
for m in agents[:Cons]
    define_common_parameters!(m,mdict[m],data,ts,agents,scenario_overview_row)                                  # Parameters common to all agents
    define_consumer_parameters!(mdict[m],merge(data["General"],data["Consumers"][m]),ts)                        # Consumers
end
# Generator models
for m in agents[:Gen]
    define_common_parameters!(m,mdict[m],data,ts,agents,scenario_overview_row)                                  # Parameters common to all agents
    define_generator_parameters!(mdict[m],merge(data["General"],data["Generators"][m]),ts)                      # Generators
end

# Calculate number of agents in each market
EOM["nAgents"] = length(agents[:eom])

println("Inititate model, sets and parameters: done")
println("   ")

## 4. Build models
for m in agents[:Cons]
    build_consumer_agent!(mdict[m])
end
for m in agents[:Gen]
    build_generator_agent!(mdict[m])
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
define_results!(merge(data["General"],data["ADMM"]),results,ADMM,agents)           # initialize structure of results, only those that will be stored in each iteration
ADMM!(results,ADMM,EOM,mdict,agents,scenario_overview_row,data,TO)                 # calculate equilibrium 
ADMM["walltime"] =  TimerOutputs.tottime(TO)*10^-9/60                              # wall time 

println(string("Done!"))
println(string("        "))
println(string("Required iterations: ",ADMM["n_iter"]))
println(string("        "))
println(string("RP EOM: ",  ADMM["Residuals"]["Primal"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
println(string("RD EOM: ",  ADMM["Residuals"]["Dual"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
println(string("        "))

## 6. Postprocessing and save results 
if sens_number >= 2
save_results(mdict,EOM,ADMM,results,data,agents,scenario_overview_row,sensitivity_overview[sens_number-1,:remarks]) 
# @save joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_",sensitivity_overview[sens_number-1,:remarks]))
else
save_results(mdict,EOM,ADMM,results,data,agents,scenario_overview_row,"ref") 
# @save joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_ref"))
end

println("Postprocessing & save results: done")
println("   ")

# end # end loop over sensititivity
# end # end for loop over scenarios

println(string("##############################################################################################"))