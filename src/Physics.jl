module Physics

using ..Constants   
using LinearAlgebra 

export newtonian_2D_model!, schwarzschild_model!


"""
classic newtonian 2D equations of motion.
u = [x, y, vx, vy]
p = mass of the black hole (M)
"""
function newtonian_2D_model!(du, u, p, t)
    # Unpack
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    M = p
    
    # Calculate r
    r = hypot(x, y) 
    
    # Newtonian Acceleration Coefficient (F/r)
    prefactor = -(Constants.G * M) / (r^3)
    
    # Update State in-place
    du[1] = vx            # dx/dt = vx
    du[2] = vy            # dy/dt = vy
    du[3] = prefactor * x # dvx/dt = ax
    du[4] = prefactor * y # dvy/dt = ay
end



"""
Includes the GR correction term from the effective potential.
"""
function schwarzschild_2D_model!(du, u, p, t)
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    M = p
    r = hypot(x, y)
    
    # Angular Momentum per unit mass 
    h = (x * vy) - (y * vx)

    # Newtonian Forces Coefficient
    term_newton = -(Constants.G * M) / r^3

    # GR Correction Coefficient
    term_gr = -(3 * Constants.G * M * h^2) / (Constants.c^2 * r^5)

    # Total Acceleration Coefficient
    total_coeff = term_newton + term_gr

    # Update State
    du[1] = vx
    du[2] = vy
    du[3] = total_coeff * x
    du[4] = total_coeff * y
end

end