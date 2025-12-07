module Visuals

using Plots
using DifferentialEquations
using Printf
using ..Constants 

export plot_orbit, animate_orbit

function plot_orbit(sol; title="Black Hole Orbit")
    x = sol[4, :]
    y = sol[5, :]
    z = sol[6, :] 
    
    M = sol.prob.p
    Rs = (2 * Constants.G * M) / Constants.c^2

    # Force Cubic Plot Volume
    max_range = maximum(abs.(vcat(x, y, z))) * 1.1 
    cube_limits = (-max_range, max_range)

    p = plot(x, y, z, 
             label="Trajectory", 
             xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
             title=title,
             aspect_ratio=:equal, 
             xlims=cube_limits,    
             ylims=cube_limits,    
             zlims=cube_limits,    
             linewidth=2,
             linecolor=:blue,
             camera=(30, 30)) 

    # Explicit Mesh Generation
    n = 20
    u = range(0, 2π, length=n)
    v = range(0, π, length=n)

    sx = [Rs * cos(U) * sin(V) for U in u, V in v]
    sy = [Rs * sin(U) * sin(V) for U in u, V in v]
    sz = [Rs * cos(V)          for U in u, V in v]

    # Draw Wireframe
    plot!(p, sx, sy, sz, color=:black, alpha=0.1, label="")
    plot!(p, sx', sy', sz', color=:black, alpha=0.1, label="")

    # Dummy label
    plot!(p, [NaN], [NaN], [NaN], color=:black, label="Event Horizon")

    scatter!(p, [x[1]], [y[1]], [z[1]], color=:green, label="Start", markersize=4)
    scatter!(p, [x[end]], [y[end]], [z[end]], color=:red, label="End", markersize=4)

    return p
end

"""
Creates a 3D GIF animation of the orbit.
"""
function animate_orbit(sol, filename="orbit.gif"; fps=30)
    M = sol.prob.p
    Rs = (2 * Constants.G * M) / Constants.c^2
    
    step_size = max(1, floor(Int, length(sol.t) / 300)) 
    indices = 1:step_size:length(sol.t)

    x_data = sol[4, :]
    y_data = sol[5, :]
    z_data = sol[6, :]
    
    max_range = maximum(abs.(vcat(x_data, y_data, z_data))) * 1.1 
    cube_limits = (-max_range, max_range)

    n = 20
    u = range(0, 2π, length=n)
    v = range(0, π, length=n)
    
    sx = [Rs * cos(U) * sin(V) for U in u, V in v]
    sy = [Rs * sin(U) * sin(V) for U in u, V in v]
    sz = [Rs * cos(V)          for U in u, V in v]

    anim = @animate for i in indices
        x_trail = sol[4, 1:i]
        y_trail = sol[5, 1:i]
        z_trail = sol[6, 1:i]
        
        p = plot(x_trail, y_trail, z_trail, 
                 xlims=cube_limits, ylims=cube_limits, zlims=cube_limits,
                 aspect_ratio=:equal,
                 label="Path", linecolor=:blue, alpha=0.5,
                 camera=(30, 30))

        scatter!(p, 
            [cube_limits[1], cube_limits[2]], 
            [cube_limits[1], cube_limits[2]], 
            [cube_limits[1], cube_limits[2]], 
            label="", alpha=0, markersize=0, color=:black
        )

        # Draw Sphere Manual Wireframe
        plot!(p, sx, sy, sz, color=:black, alpha=0.15, label="")
        plot!(p, sx', sy', sz', color=:black, alpha=0.15, label="")
        
        # Draw current position
        scatter!(p, [sol[4, i]], [sol[5, i]], [sol[6, i]], 
                 color=:red, label="Particle", markersize=4)
        
        title!(p, "Time: $(round(sol.t[i], digits=2)) s")
    end

    gif(anim, filename, fps=fps)
    println("Animation saved to $filename")
end

end