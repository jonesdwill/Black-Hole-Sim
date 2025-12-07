module Solvers

using DifferentialEquations
using LinearAlgebra
using RecursiveArrayTools
using ..Constants
using ..Physics

export setup_problem, solve_orbit

function setup_problem(model_type::Symbol, u0, tspan, M)
    pos_0 = u0[1:3] 
    vel_0 = u0[4:6] 

    if model_type == :schwarzschild

        return DynamicalODEProblem(Physics.schwarzschild_acceleration!, Physics.velocity_law!, vel_0, pos_0, tspan, M)

    elseif model_type == :newtonian

        return DynamicalODEProblem(Physics.newtonian_acceleration!, Physics.velocity_law!, vel_0, pos_0, tspan, M)

    else
        error("Unknown model type.")
    end
end

function solve_orbit(prob; rel_tol=1e-12, dt=1e-5, kwargs...)
    
    M = prob.p
    Rs = (2 * Constants.G * M) / Constants.c^2

    function condition(u, t, integrator)
        pos = u.x[2]      
        r = norm(pos)
        return r - Rs 
    end

    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    return solve(prob, KahanLi8(), dt=dt, reltol=rel_tol, abstol=rel_tol, callback=cb)
end

end