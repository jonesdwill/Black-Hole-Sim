module BasicBlackHoleSim

# External Packages
using DifferentialEquations
using Plots
using LinearAlgebra
using RecursiveArrayTools

# Project Structure 
include("Constants.jl")
include("Physics.jl")
include("Solvers.jl")
include("Visuals.jl")

# Bring into scope if needed
using .Constants
using .Physics
using .Solvers
using .Visuals

# User-end
export Constants                   # Access G, c
export NewtonianPhysics            # Access equations
export simulate_orbit              # Main
export plot_orbit, animate_orbit   # Visualisation


"""
High-level wrapper to run a simulation.
"""
function simulate_orbit(model_type, u0, tspan, M)
    prob = Solvers.setup_problem(model_type, u0, tspan, M)
    return Solvers.solve_orbit(prob)
end

export simulate_orbit

end 