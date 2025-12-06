module Solvers

using DifferentialEquations
using LinearAlgebra

using ..Constants
using ..Physics

export setup_problem, solve_orbit

"""
Creates an ODE for the specific physics model.
"""
function setup_problem(model_type::Symbol, u0, tspan, M)
    
    if model_type == :newtonian
        func = Physics.newtonian_2D_model!
    elseif model_type == :schwarzschild
        func = Physics.schwarzschild_2D_model!
    else
        error("Unknown model type. Use :newtonian or :schwarzschild")
    end

    return ODEProblem(func, u0, tspan, M)
end

"""
Solves ODE with high precision.
Includes a callback to stop if the particle hits Event Horizon.
"""
function solve_orbit(prob; rel_tol=1e-12)
    
    # safety check
    M = prob.p
    Rs = (2 * G * M) / c^2

    function condition(u, t, integrator)
        r = norm(u[1:2])
        return r - Rs # When this hits 0, the event triggers
    end

    # stop 
    affect!(integrator) = terminate!(integrator)
    
    # finds when r == Rs
    cb = ContinuousCallback(condition, affect!)

    # Vern9 is a highly accurate 9th-order Runge-Kutta solver.
    return solve(prob, Vern9(), reltol=rel_tol, abstol=rel_tol, callback=cb)
end

end