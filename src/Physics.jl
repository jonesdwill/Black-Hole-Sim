module Physics

using ..Constants 
using LinearAlgebra
using ForwardDiff

# EXPORTS: We now export the split functions
export velocity_law!, newtonian_acceleration!, schwarzschild_acceleration!, kerr_acceleration!, kerr_geodesic!

# ===========================================
#       POST-NEWTONIAN APPROXIMATIONS
# ===========================================

"""
VELOCITY LAW
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

    local G = Constants.G
    local M = p

    x, y, z = q[1], q[2], q[3]
    r = norm(q)

    prefactor = -(G * M) / (r^3)
    
    # Update Acceleration (dv)
    dv[1] = prefactor * x
    dv[2] = prefactor * y
    dv[3] = prefactor * z
end

"""
SCHWARZSCHILD ACCELERATION
"""
function schwarzschild_acceleration!(dv, v, q, p, t)
    
    local G = Constants.G       # gravity constant
    local M = p                 # mass
    local c = Constants.c       # speed of light

    x, y, z = q[1], q[2], q[3]  # pos
    r_sq = x^2 + y^2 + z^2      # radius distance norn 
    r = sqrt(r_sq)              # radius 
    
    h_vec = cross(q, v) 
    h_sq = dot(h_vec, h_vec) 

    # Coefficients
    term_newton = -(G * M) / (r * r_sq)                           # Newtonian term: -GM/r^3
    term_gr = -(3 * G * M * h_sq) / (c^2 * r_sq^2 * r)  # General Relativity Term: -3GMh^2/(c^2 r^5)

    # Total force applied 
    total_coeff = term_newton + term_gr

    # Update Acceleration
    dv[1] = total_coeff * x
    dv[2] = total_coeff * y
    dv[3] = total_coeff * z 
end


"""
KERR ACCELERATION 
Includes Schwarzschild precession and Lense-Thirring frame-dragging.
"""
function kerr_acceleration!(dv, v, q, p, t)
    
    local G = Constants.G       # gravity constant
    local M = p                 # mass
    local c = Constants.c       # speed of light

    x, y, z = q[1], q[2], q[3]
    vx, vy, vz = v[1], v[2], v[3]
    M, a_star = p 
    r_sq = x^2 + y^2 + z^2
    r = sqrt(r_sq)

    # --- Schwarzschild part (central force) ---
    h_vec = cross(q, v) 
    h_sq = dot(h_vec, h_vec) 

    term_newton = -(G * M) / (r * r_sq) # -GM/r^3
    term_gr = -(3 * G * M * h_sq) / (c^2 * r_sq^2 * r) # -3GMh^2/(c^2 r^5)

    central_coeff = term_newton + term_gr

    ax_central = central_coeff * x
    ay_central = central_coeff * y
    az_central = central_coeff * z

    # --- Lense-Thirring part (frame-dragging, non-central force) ---
    if a_star != 0.0

        # black hole's angular momentum
        J_mag = a_star * G * M^2 / c
        
        # Prefactor Lense-Thirring acceleration
        prefactor_lt = (2 * G * J_mag) / (c^2 * r * r_sq)
        
        # Lense-Thirring acceleration components, in newtonian approximation
        ax_lt = prefactor_lt * ( (3 * z / r_sq) * h_vec[1] + vy )
        ay_lt = prefactor_lt * ( (3 * z / r_sq) * h_vec[2] - vx )
        az_lt = prefactor_lt * ( (3 * z / r_sq) * h_vec[3] )

        # Total acceleration
        dv[1] = ax_central + ax_lt
        dv[2] = ay_central + ay_lt
        dv[3] = az_central + az_lt

    else
        # If a_star is 0, just use Schwarzschild
        dv[1] = ax_central
        dv[2] = ay_central
        dv[3] = az_central
    end
end

# ===========================================
#           GEOMETRIC APPROXIMATION
# ===========================================

# #  Hard coded geodesic Kerr equations here: nothing spectacular

"""
KERR ACCELERATION 
Includes Schwarzschild precession and Lense-Thirring frame-dragging.
"""
function kerr_geodesic_acceleration!(du, u, p, λ)
    M, a = p
    
    # Coordinates
    t, r, theta, phi = u[1], u[2], u[3], u[4]
    # Momenta (covariant components p_mu)
    pt, pr, ptheta, pphi = u[5], u[6], u[7], u[8]

    # --- Metric Inverse Components (g^muν) ---
    # Sigma = r^2 + a^2 cos^2(theta)
    # Delta = r^2 - 2Mr + a^2
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sigma = r^2 + a^2 * cos_theta^2
    delta = r^2 - 2*M*r + a^2

    # Contravariant Metric components g^muν
    inv_g_tt = -((r^2 + a^2)^2 - delta * a^2 * sin_theta^2) / (delta * sigma)
    inv_g_rr = delta / sigma
    inv_g_thetatheta = 1.0 / sigma
    inv_g_phiphi = (delta - a^2 * sin_theta^2) / (delta * sigma * sin_theta^2)
    inv_g_tphi = -(2 * M * r * a) / (delta * sigma)

    # --- Hamilton's Equations ---
    # 1Velocity equations: dx^mu / dλ = g^muν p_ν
    
    # dt/dλ
    du[1] = inv_g_tt * pt + inv_g_tphi * pphi
    # dr/dλ
    du[2] = inv_g_rr * pr
    # dtheta/dλ
    du[3] = inv_g_thetatheta * ptheta
    # dphi/dλ
    du[4] = inv_g_tphi * pt + inv_g_phiphi * pphi

    # Force equations: dp_mu / dλ = -1/2 * (∂g^αβ / ∂x^mu) p_α p_β

    # --- Derivative wrt r ---
    # ∂sigma/∂r = 2r, ∂delta/∂r = 2r - 2M
    dsigma_dr = 2*r
    ddelta_dr = 2*r - 2*M
    
    # Define H for Automatic Differentiation 
    H_func = (coords) -> begin
        _r, _theta = coords[1], coords[2]
        _stheta = sin(_theta); _ctheta = cos(_theta)
        _sigma = _r^2 + a^2 * _ctheta^2
        _delta = _r^2 - 2*M*_r + a^2
        
        _gtt = -((_r^2 + a^2)^2 - _delta * a^2 * _stheta^2) / (_delta * _sigma)
        _grr = _delta / _sigma
        _gthth = 1.0 / _sigma
        _gphph = (_delta - a^2 * _stheta^2) / (_delta * _sigma * _stheta^2)
        _gtph = -(2 * M * _r * a) / (_delta * _sigma)
        
        return 0.5 * (_gtt*pt^2 + _grr*pr^2 + _gthth*ptheta^2 + _gphph*pphi^2 + 2*_gtph*pt*pphi)
    end

    # Compute gradients for force terms
    grads = ForwardDiff.gradient(H_func, [r, theta])

    du[5] = 0.0          # dp_t / dλ (Constant of motion E)
    du[6] = -grads[1]    # dp_r / dλ = -dH/dr
    du[7] = -grads[2]    # dp_theta / dλ = -dH/dtheta
    du[8] = 0.0          # dp_phi / dλ (Constant of motion Lz)
end

end