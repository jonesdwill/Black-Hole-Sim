using Test
using BasicBlackHoleSim
using BasicBlackHoleSim.Constants
using BasicBlackHoleSim.Utils
using BasicBlackHoleSim.Solvers
using LinearAlgebra
using DifferentialEquations 

# --- Test Helper Functions ---

"""
Calculates the specific orbital energy (energy per unit mass) for a Newtonian orbit.
E/m = 1/2 v^2 - G M / r
"""
function calculate_newtonian_energy(u, M)
    q = u[1:3]
    v = u[4:6]
    
    r = norm(q)
    v_sq = dot(v, v)
    
    kinetic_energy = 0.5 * v_sq
    potential_energy = -(G * M) / r
    
    return kinetic_energy + potential_energy
end

"""
Calculates the Hamiltonian for a particle in Kerr spacetime.
H = 1/2 g^μν p_μ p_ν
For a massive particle (m=1), this should be conserved and equal to -0.5.
"""
function calculate_kerr_hamiltonian(u, p)
    M, a = p
    
    # Coordinates and Momenta
    r, θ = u[2], u[3]
    pt, pr, pθ, pϕ = u[5], u[6], u[7], u[8]

    # Metric Inverse Components (g^μν)
    sin_θ, cos_θ = sincos(θ)
    Σ = r^2 + a^2 * cos_θ^2
    Δ = r^2 - 2*M*r + a^2

    # Avoid division by zero if at singularity or on horizon with certain coordinates
    if Σ == 0 || Δ == 0 || sin_θ == 0
        return NaN
    end

    inv_g_tt = -((r^2 + a^2)^2 - Δ * a^2 * sin_θ^2) / (Δ * Σ)
    inv_g_rr = Δ / Σ
    inv_g_θθ = 1.0 / Σ
    inv_g_ϕϕ = (Δ - a^2 * sin_θ^2) / (Δ * Σ * sin_θ^2)
    inv_g_tϕ = -(2 * M * r * a) / (Δ * Σ)

    # Hamiltonian H = 1/2 * g^μν p_μ p_ν
    return 0.5 * (inv_g_tt*pt^2 + inv_g_rr*pr^2 + inv_g_θθ*pθ^2 + inv_g_ϕϕ*pϕ^2 + 2*inv_g_tϕ*pt*pϕ)
end

@testset "BasicBlackHoleSim.jl" begin

    @testset "Physics Helpers" begin
        M = Constants.M # Use dimensionless mass (1.0)
        
        @testset "get_black_hole_parameters" begin
            expected_Rs = (2 * G * M) / c^2
            params = get_black_hole_parameters(M)
            @test params.rs ≈ expected_Rs
        end
        
        @testset "circular_velocity" begin
            r0 = 5.0e7
            expected_v = sqrt(G * M / r0)
            @test circular_velocity(M, r0) ≈ expected_v
        end

        @testset "get_initial_photon_state_scattering" begin
            M = 1.0
            a = 0.9
            r0 = 50.0
            b = 5.0 # impact parameter
            u0_photon = get_initial_photon_state_scattering(r0, b, M, a)
            
            # For a photon, the Hamiltonian must be 0
            H_photon = calculate_kerr_hamiltonian(u0_photon, (M, a))
            @test H_photon ≈ 0.0 atol=1e-9
        end
    end

    @testset "Post-Newtonian Solvers" begin
        M = Constants.M # Use dimensionless mass (1.0)
        r0 = 20.0 # Use a larger radius for stability in PN tests
        v0 = circular_velocity(M, r0)
        tspan = (0.0, 100.0) # Longer simulation time for better test of conservation

        u0 = [r0, 0.0, 0.0, 0.0, v0, 0.0] 

        @testset "Newtonian Solver" begin
            sol = simulate_orbit(:newtonian, u0, tspan, M)
            @test Symbol(sol.retcode) == :Success
            
            # Test for energy conservation
            E_initial = calculate_newtonian_energy(sol.u[1], M)
            E_final = calculate_newtonian_energy(sol.u[end], M)
            @test E_final ≈ E_initial rtol=1e-6 # Check conservation to a reasonable tolerance
        end

        @testset "Schwarzschild Solver" begin
            sol = simulate_orbit(:schwarzschild, u0, tspan, M)
            @test Symbol(sol.retcode) == :Success
            @test sol.u[end] != u0
        end

        @testset "Kerr Solver" begin
            a_star = 0.98
            p = (M, a_star)
            sol = simulate_orbit(:kerr_acceleration, u0, tspan, p)
            @test Symbol(sol.retcode) == :Success
            @test sol.u[end] != u0
        end
    end

    @testset "Geodesic Kerr Solver" begin
        M = 1.0
        a_star = 0.9
        kerr_params = (M, a_star)
        # In geometric units with M=1, a_geom = a_star
        a_geom = a_star 

        r0 = 10.0
        theta0 = π/2
        phi0 = 0.0
        tspan = (0.0, 500.0)

        # Use utils to get initial conditions for a circular orbit
        ut_circ, uphi_circ = calculate_circular_orbit_properties(r0, M, a_geom)
        
        # Geodesic state vector [t, r, θ, ϕ, uᵗ, uʳ, uᶿ, uᵠ]
        u_geo_initial = [0.0, r0, theta0, phi0, ut_circ, 0.0, 0.0, uphi_circ]

        # Convert to Hamiltonian state vector [t, r, θ, ϕ, pₜ, pᵣ, pₜ, pᵩ]
        u_ham_initial = get_initial_hamiltonian_state(u_geo_initial, M, a_geom)

        # Setup and solve
        prob = setup_problem(:kerr_geodesic_acceleration, u_ham_initial, tspan, kerr_params)
        sol = solve_orbit(prob, reltol=1e-10, abstol=1e-10) # Use tighter tolerances for conservation tests

        @test Symbol(sol.retcode) == :Success

        # Test conservation of energy (p_t) and angular momentum (p_phi)
        pt_initial = u_ham_initial[5]
        p_phi_initial = u_ham_initial[8]
        @test sol.u[end][5] ≈ pt_initial rtol=1e-7
        @test sol.u[end][8] ≈ p_phi_initial rtol=1e-7

        # Test conservation of Hamiltonian (related to rest mass)
        H_initial = calculate_kerr_hamiltonian(sol.u[1], kerr_params)
        H_final = calculate_kerr_hamiltonian(sol.u[end], kerr_params)

        # For a massive particle with m=1, H = -m^2/2 = -0.5
        @test H_initial ≈ -0.5 rtol=1e-7
        @test H_final ≈ H_initial rtol=1e-7
    end
    
end