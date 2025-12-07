module Physics

using ..Constants 
using LinearAlgebra 

# EXPORTS: We now export the split functions
export velocity_law!, newtonian_acceleration!, schwarzschild_acceleration!

"""
THE VELOCITY LAW
"""
function velocity_law!(dq, v, q, p, t)
    dq[1] = v[1]
    dq[2] = v[2]
    dq[3] = v[3]
end

"""
NEWTONIAN ACCELERATION
"""
function newtonian_acceleration!(dv, v, q, p, t)

    x, y, z = q[1], q[2], q[3]
    M = p
    r = norm(q)

    prefactor = -(Constants.G * M) / (r^3)
    
    # Update Acceleration (dv)
    dv[1] = prefactor * x
    dv[2] = prefactor * y
    dv[3] = prefactor * z
end

"""
SCHWARZSCHILD ACCELERATION
"""
function schwarzschild_acceleration!(dv, v, q, p, t)
    x, y, z = q[1], q[2], q[3]
    vx, vy = v[1], v[2]
    M = p
    r = norm(q)
    
    h_vec = cross(q, v) 
    h_sq = dot(h_vec, h_vec) 

    # Coefficients
    term_newton = -(Constants.G * M) / r^3
    term_gr = -(3 * Constants.G * M * h_sq) / (Constants.c^2 * r^5)

    total_coeff = term_newton + term_gr

    # Update Acceleration
    dv[1] = total_coeff * x
    dv[2] = total_coeff * y
    dv[3] = total_coeff * z 
end

end