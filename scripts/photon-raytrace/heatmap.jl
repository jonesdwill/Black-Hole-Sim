# projectdir(args...) = joinpath(@__DIR__, "..", "..", args...)
# import Pkg; Pkg.activate(projectdir())

# using BasicBlackHoleSim
# using Plots
# using BasicBlackHoleSim.Solvers: setup_problem, solve_orbit
# using BasicBlackHoleSim.Utils: get_black_hole_parameters, get_initial_photon_state_scattering
# using Printf
# using Base.Threads # To speed up calculation

# # ==========================================
# # 1. SETUP
# # ==========================================
# M = 1.0
# a_star = 0.99 
# solver_params = (M, a_star) 
# a = a_star * M
# bh_params = get_black_hole_parameters(solver_params)

# # Simulation Settings
# # Heatmaps need MANY more particles to look smooth, 
# # but because we aren't drawing lines, we can handle thousands easily.
# n_particles = 2000      
# max_life = 60.0        
# anim_duration = 30.0   
# emission_rate = 5.0     # How often (in time units) a new wave emits

# # View / Resolution
# x_view = (-12, 12)
# y_view = (-8, 8)
# heatmap_bins = 200      # Resolution of the heatmap grid (200x200 pixels)

# # ==========================================
# # 2. PARALLEL PRE-COMPUTATION
# # ==========================================
# println("1/3: Calculating $(n_particles) geodesics on $(Threads.nthreads()) threads...")

# # We store the solutions directly to sample them later
# solutions = Vector{Any}(undef, n_particles)
# b_values = range(-10.0, 10.0, length=n_particles)

# # Use multi-threading to speed up the math
# Threads.@threads for i in 1:n_particles
#     b = b_values[i]
    
#     # Skip center singularity to avoid NaNs
#     if abs(b) < 0.1 
#         solutions[i] = nothing
#         continue 
#     end

#     u0 = get_initial_photon_state_scattering(40.0, b, M, a)
#     prob = setup_problem(:kerr_geodesic_acceleration, u0, (0.0, max_life * 1.5), solver_params)
    
#     # Lower tolerance slightly for speed since heatmaps blur details anyway
#     solutions[i] = solve_orbit(prob, reltol=1e-5, abstol=1e-5)
# end

# # Filter out failed solves
# valid_sols = filter(!isnothing, solutions)
# println("     Computed $(length(valid_sols)) valid orbits.")

# # ==========================================
# # 3. ANIMATION
# # ==========================================
# println("2/3: Rendering heatmap frames...")

# # Pre-calculate the Black Hole Shadow for drawing the circle
# u_circ = range(0, 2π, length=100)
# hx = bh_params.rh .* cos.(u_circ)
# hy = bh_params.rh .* sin.(u_circ)

# anim = @animate for t_global in range(0, anim_duration, step=0.2)
    
#     # Collect all X and Y positions for this specific moment in time
#     # We simulate "waves" of particles being emitted repeatedly
#     xs = Float64[]
#     ys = Float64[]
    
#     num_waves = floor(Int, max_life / emission_rate)
    
#     for i in 0:num_waves
#         # Calculate the "local time" for this wave
#         # This creates the continuous flow effect
#         local_t = (t_global % emission_rate) + (i * emission_rate)
        
#         # Optimization: Don't sample if outside bounds
#         if local_t > max_life || local_t < 0
#             continue
#         end

#         # Sample every pre-calculated orbit at this time
#         for sol in valid_sols
#             # Check if particle is still alive (hasn't hit singularity/max time)
#             if local_t > sol.t[end]
#                 continue
#             end

#             state = sol(local_t)
#             r, th, ph = state[2], state[3], state[4]

#             # Filter Horizon hits (don't draw points inside the black hole)
#             if r < bh_params.rh * 1.05
#                 continue
#             end

#             # Convert Boyer-Lindquist to Cartesian for plotting
#             xx = sqrt(r^2 + a^2) * sin(th) * cos(ph)
#             yy = sqrt(r^2 + a^2) * sin(th) * sin(ph)

#             # Viewport culling (don't process points outside the camera view)
#             if x_view[1] < xx < x_view[2] && y_view[1] < yy < y_view[2]
#                 push!(xs, xx)
#                 push!(ys, yy)
#             end
#         end
#     end

#     # --- THE PLOT ---
#     # histogram2d does the heavy lifting: bins points into a grid
#     histogram2d(
#         xs, ys,
#         bins = heatmap_bins,
#         xlims = x_view,
#         ylims = y_view,
#         c = :inferno,        # Fire colormap
#         bg = :black,
#         clims = (0, 15),     # Clamp brightness so the hole edge glows but doesn't wash out
#         legend = false,
#         aspect_ratio = :equal,
#         framestyle = :none,
#         margin = -5Plots.mm
#     )

#     # Draw the black hole shadow on top
#     plot!(hx, hy, seriestype=:shape, c=:black, lw=0)
    
# end

# # ==========================================
# # 4. SAVE
# # ==========================================
# println("3/3: Saving GIF...")
# output = projectdir("scripts", "kerr_heatmap.gif")
# gif(anim, output, fps = 20)
# println("Saved: $output")

projectdir(args...) = joinpath(@__DIR__, "..", "..", args...)
import Pkg; Pkg.activate(projectdir())

using BasicBlackHoleSim
using Plots
using BasicBlackHoleSim.Solvers: setup_problem, solve_orbit
using BasicBlackHoleSim.Utils: get_black_hole_parameters, get_initial_photon_state_scattering
using Random
using Base.Threads

# ==========================================
# 1. HIGH-QUALITY SETUP
# ==========================================
M = 1.0
a_star = 0.99 
solver_params = (M, a_star) 
a = a_star * M
bh_params = get_black_hole_parameters(solver_params)

# Rendering Settings
n_particles = 10_000    # Huge increase for density
max_life = 60.0        
anim_duration = 30.0   
emission_rate = 5.0     
heatmap_bins = 400      # Higher resolution grid (400x400)

# "Smear" settings to fill gaps
samples_per_streak = 5  # How many points to draw per particle per frame
streak_length = 0.4     # Length of the trail (should match animation step size)

x_view = (-12, 12)
y_view = (-8, 8)

# ==========================================
# 2. PARALLEL CALCULATION
# ==========================================
println("1/3: Calculating $(n_particles) geodesics (High Quality)...")

solutions = Vector{Any}(undef, n_particles)

# Add randomness to 'b' to break the artificial vertical lines
# This fills the "spatial gaps" between the streams
rng = range(-10.0, 10.0, length=n_particles)
b_values = [b + (rand() * 0.05) for b in rng] 

Threads.@threads for i in 1:n_particles
    b = b_values[i]
    if abs(b) < 0.1 
        solutions[i] = nothing
        continue 
    end

    u0 = get_initial_photon_state_scattering(40.0, b, M, a)
    prob = setup_problem(:kerr_geodesic_acceleration, u0, (0.0, max_life * 1.5), solver_params)
    solutions[i] = solve_orbit(prob, reltol=1e-5, abstol=1e-5)
end

valid_sols = filter(!isnothing, solutions)
println("     Computed $(length(valid_sols)) valid orbits.")

# ==========================================
# 3. ANIMATION WITH STREAKING
# ==========================================
println("2/3: Rendering smooth frames...")

u_circ = range(0, 2π, length=100)
hx = bh_params.rh .* cos.(u_circ)
hy = bh_params.rh .* sin.(u_circ)

# We step time by 'streak_length' so the streaks connect perfectly
anim = @animate for t_global in range(0, anim_duration, step=0.2)
    
    xs = Float64[]
    ys = Float64[]
    
    num_waves = floor(Int, max_life / emission_rate)
    
    # Pre-allocate somewhat to avoid resizing constantly
    sizehint!(xs, 100_000)
    sizehint!(ys, 100_000)

    for i in 0:num_waves
        base_t = (t_global % emission_rate) + (i * emission_rate)
        
        if base_t > max_life || base_t < 0 continue end

        Threads.@threads for sol in valid_sols
            # --- STREAK LOGIC ---
            # Instead of sampling 1 point, sample 'samples_per_streak' points 
            # backward from the current time. This draws a line, filling the gap.
            for k in 1:samples_per_streak
                # Calculate a sub-time slightly in the past
                sub_t = base_t - ((k-1) / samples_per_streak) * streak_length
                
                if sub_t > sol.t[end] || sub_t < 0 continue end

                state = sol(sub_t)
                r, th, ph = state[2], state[3], state[4]

                if r < bh_params.rh * 1.05 continue end

                xx = sqrt(r^2 + a^2) * sin(th) * cos(ph)
                yy = sqrt(r^2 + a^2) * sin(th) * sin(ph)

                if x_view[1] < xx < x_view[2] && y_view[1] < yy < y_view[2]
                    # We use a lock because push! is not thread-safe
                    # (Or we can just do this loop single-threaded if it errors, usually it's fast enough)
                    lock(ReentrantLock()) do
                        push!(xs, xx)
                        push!(ys, yy)
                    end
                end
            end
        end
    end

    histogram2d(
        xs, ys,
        bins = heatmap_bins,
        xlims = x_view,
        ylims = y_view,
        c = :inferno,
        bg = :black,
        clims = (0, 30), # Higher limit because we have more overlapping points now
        legend = false,
        aspect_ratio = :equal,
        framestyle = :none,
        margin = -5Plots.mm
    )

    plot!(hx, hy, seriestype=:shape, c=:black, lw=0)
    
end

# ==========================================
# 4. SAVE
# ==========================================
println("3/3: Saving High-Res GIF...")
output = projectdir("scripts", "kerr_heatmap_smooth.gif")
gif(anim, output, fps = 24)
println("Saved: $output")