projectdir(args...) = joinpath(@__DIR__, "..", "..", args...)
import Pkg; Pkg.activate(projectdir())

using Plots
using Serialization
using Printf

# ==========================================
# 1. ANIMATION SETTINGS
# ==========================================
num_frames = 60
fps = 15

# Visuals
base_exposure = 3.5
rotation_speed = 6.0   # How fast it spins
swirl_density = 5.0    # Lower = Smoother, fatter spiral arms (Less glitchy)

function get_pixel_color(r, redshift, phi, time_val)
    if r < 0 
        return RGB(0.0, 0.0, 0.0) 
    end

    # 1. Doppler Beaming (Static)
    beaming = clamp(redshift^4, 0.05, 50.0) 
    
    # 2. SPINNING TEXTURE
    # We add 'time_val' to 'phi'. Since phi is the angle around the circle,
    # adding to it rotates the pattern!
    noise = sin(phi * swirl_density - time_val * rotation_speed + r)
    
    # Map noise to brightness
    texture_variation = 1.0 + 0.4 * noise
    
    intensity = base_exposure * beaming * texture_variation

    # 3. Fire Palette
    if redshift > 1.1
        return RGB(min(1.0, intensity*0.8), min(1.0, intensity*0.9), min(1.0, intensity))
    elseif redshift > 0.8
        return RGB(min(1.0, intensity), min(1.0, intensity*0.6), 0.0)
    else
        return RGB(min(1.0, intensity*0.5), 0.0, 0.0)
    end
end

# ==========================================
# 2. RENDER LOOP
# ==========================================
data_file = projectdir("scripts/backwards-raytrace", "static_blackhole_data.jls")

if !isfile(data_file)
    println("Run the Compute script first!")
    exit()
end

println("Loading Data...")
raw_data = deserialize(data_file)
height, width = size(raw_data)

println("Rendering Spinning Animation...")

anim = @animate for frame_idx in 1:num_frames
    
    # Calculate Time (0.0 to 1.0)
    progress = (frame_idx - 1) / num_frames
    current_time = progress * 2π # Complete loop
    
    image_grid = fill(RGB(0.0, 0.0, 0.0), height, width)

    for j in 1:height
        for i in 1:width
            r_val, g_val, phi_val = raw_data[j, i]
            image_grid[j, i] = get_pixel_color(r_val, g_val, phi_val, current_time)
        end
    end
    
    plot(image_grid, axis=nothing, border=:none, size=(width*2, height*2), aspect_ratio=:equal)
    print("\rCompositing Frame $frame_idx / $num_frames")
end

output_path = projectdir("scripts/backwards-raytrace", "blackhole_spin.gif")
gif(anim, output_path, fps=fps)
println("\nSaved to: $output_path")