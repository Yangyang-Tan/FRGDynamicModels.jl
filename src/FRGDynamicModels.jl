module FRGDynamicModels
using Sundials,ProgressLogging
using DifferentialEquations
using SpecialFunctions
using Dierckx
using Plots
include("flow.jl")
include("solver.jl")
export solver_modelA
export plot
end
