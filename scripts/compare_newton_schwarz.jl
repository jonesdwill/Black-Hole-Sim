using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 
using BasicBlackHoleSim
using Plots

# --- CONFIGURATION ---
M = 2.0e30  
tspan = (0.0, 0.2) 

u0 = [100.0e3, 0.0, 0.0, 3.6e7] 

println("1. Running Newtonian Simulation...")
sol_newton = BasicBlackHoleSim.Solvers.solve_orbit(
    BasicBlackHoleSim.Solvers.setup_problem(:newtonian, u0, tspan, M),
    dt=1e-5
)

println("2. Running Schwarzschild Simulation...")
sol_gr = BasicBlackHoleSim.Solvers.solve_orbit(
    BasicBlackHoleSim.Solvers.setup_problem(:schwarzschild, u0, tspan, M),
    dt=1e-5
)

println("3. Generating Comparison Plot...")

# Plot Newtonian 
p = plot(sol_newton[3, :], sol_newton[4, :], 
         label="Newtonian (Standard)", color=:green, linewidth=1,
         aspect_ratio=:equal, title="Gravity Comparison: Newton vs Einstein")

# Plot Schwarzschild 
plot!(p, sol_gr[3, :], sol_gr[4, :], 
      label="Schwarzschild (GR)", color=:blue, linewidth=1, alpha=0.7)

# Add Black Hole
Rs = (2 * 6.6743e-11 * M) / (2.9979e8^2)
theta = range(0, 2π, length=100)
plot!(p, Rs .* cos.(theta), Rs .* sin.(theta), 
      seriestype=[:shape], c=:black, label="Event Horizon")

display(p)
savefig(p, joinpath(@__DIR__, "comparison_result.png"))
println("   -> Saved comparison to comparison_result.png")