projectdir(args...) = joinpath(@__DIR__, "..", "..", args...)
import Pkg; Pkg.activate(projectdir())

using BasicBlackHoleSim
using DifferentialEquations
using LinearAlgebra
using Base.Threads
using Serialization
using Printf

# Access internal modules
using BasicBlackHoleSim.Constants
using BasicBlackHoleSim.Physics
using BasicBlackHoleSim.Utils

# ==========================================
# 1. CONFIG (One Static Frame)
# ==========================================
width = 400
height = 225

# Fixed Camera Angle (75 degrees looks cool)
r_cam = 1000.0
theta_cam = 75.0 * (π/180) 
FOV = 25.0

M = 1.0
a_star = 0.99
a = a_star * M
rh = M + sqrt(M^2 - a^2)
r_isco = 1.0 + sqrt(1.0 - a_star^2)
disk_outer_edge = 20.0

# ==========================================
# 2. PHYSICS ENGINE
# ==========================================
function get_redshift(r, p_t, p_phi, M, a)
    Omega = 1.0 / (r^1.5 + a)
    g_tt = -(1 - 2*M/r)
    g_tphi = -2*M*a/r
    g_phiphi = (r^2 + a^2 + 2*M*a^2/r)
    
    norm_sq = g_tt + 2*Omega*g_tphi + Omega^2*g_phiphi
    if norm_sq >= 0 return 0.0 end
    u_t = 1.0 / sqrt(-norm_sq)

    E_inf = -p_t
    E_emit = -(p_t * u_t + p_phi * (u_t * Omega))
    return E_inf / E_emit
end

println("Computing Static Frame ($width x $height)...")

raw_data = fill((-1.0f0, 0.0f0, 0.0f0), height, width)

alphas = range(-FOV/2, FOV/2, length=width)
betas = range(-FOV*(height/width)/2, FOV*(height/width)/2, length=height)
params = (M, a)

horizon_cb = ContinuousCallback((u,t,i)->u[2]-(rh*1.01), terminate!)
disk_cb = ContinuousCallback((u,t,i)->cos(u[3]), terminate!)
cb_set = CallbackSet(horizon_cb, disk_cb)

counter = Threads.Atomic{Int}(0)
total_pixels = width * height

Threads.@threads for j in 1:height
    for i in 1:width
        alpha = alphas[i]
        beta = betas[j]

        u0_vec = get_initial_photon_state_celestial(alpha, beta, r_cam, theta_cam, M, a)
        if any(isnan, u0_vec) continue end

        # High Quality Settings (reltol=1e-4) to prevent "static noise"
        prob = ODEProblem(Physics.kerr_geodesic_acceleration!, u0_vec, (0.0, 3000.0), params)
        sol = solve(prob, Tsit5(), reltol=1e-4, abstol=1e-4, callback=cb_set, save_everystep=false, dtmax=20.0)

        final_state = sol.u[end]
        r_final = final_state[2]
        theta_final = final_state[3]
        phi_final = final_state[4] 
        
        hit_disk = abs(cos(theta_final)) < 0.05 && r_final > rh * 1.05
        
        if hit_disk && r_final > r_isco && r_final < disk_outer_edge
            p_t = final_state[5]
            p_phi = final_state[8]
            g = get_redshift(r_final, p_t, p_phi, M, a)
            
            # Save R, G, and PHI
            raw_data[j, i] = (Float32(r_final), Float32(g), Float32(phi_final))
        end

        c = Threads.atomic_add!(counter, 1)
        if c % 5000 == 0
            print("\rProgress: $(round(c/total_pixels*100, digits=1))% ")
        end
    end
end

output_file = projectdir("scripts/backwards-raytrace", "static_blackhole_data.jls")
serialize(output_file, raw_data)
println("\nDone! Saved to $output_file")