# Kerr Black Hole Accretion Disk Renderer

This project is a physically-based renderer written in Julia to simulate and visualize the appearance of an accretion disk around a rotating (Kerr) black hole. It uses numerical geodesic integration to trace light paths and includes relativistic effects like gravitational lensing, Doppler beaming, and gravitational redshift.

![Accretion Disk Animation](scripts/8-full-render/black_hole_1080p.gif)

## Features

*   **Kerr Spacetime:** Simulates the geometry around a rotating black hole using the Kerr metric.
*   **Geodesic Ray Tracing:** Traces the paths of photons from a virtual camera back in time to determine what is seen at each pixel.
*   **Physically-Motivated Accretion Disk:** Models a thin accretion disk orbiting in the equatorial plane, bounded by the Innermost Stable Circular Orbit (ISCO).
*   **Relativistic Effects Visualized:**
    *   **Gravitational Lensing:** The bending of light from the disk as it passes near the black hole.
    *   **Doppler Beaming:** The side of the disk moving towards the camera appears significantly brighter and blueshifted.
    *   **Gravitational Redshift:** Light from deeper in the gravitational well is redshifted.
    *   **Black Hole Shadow:** The dark region where light rays are captured by the event horizon.
*   **Advanced Procedural Shader:**
    *   Multi-layered noise functions create a turbulent, gaseous appearance.
    *   Spiral density waves propagate through the disk.
    *   Color is determined by the calculated redshift factor `g`.
*   **High-Quality Rendering Pipeline:**
    *   A two-step process separates the expensive physics calculations from the faster rendering step.
    *   Post-processing effects like bloom (glow) and tone mapping are applied for a cinematic look.
    *   Outputs to high-resolution images or animated GIFs.

## How It Works

The simulation is split into two main parts, controlled by flags in the render script:

1.  **Physics Computation (`RUN_PHYSICS = true`):** For each pixel on the screen, a photon's path (a null geodesic) is traced backwards in time. The simulation solves the geodesic equations until the photon either hits the accretion disk, falls into the event horizon, or escapes. If it hits the disk, the intersection point and the redshift factor (`g`) are saved to a `.jld2` data file. This step is computationally intensive and only needs to be run once for a given black hole and camera configuration.

2.  **Animation Rendering (`RUN_RENDER = true`):** This step loads the pre-computed data file. For each frame of the animation, it calculates the appearance of the disk by applying the procedural shader. It rotates the disk, applies turbulence and color based on the saved data, and performs post-processing. This step is much faster and can be run repeatedly to tweak visual parameters like color, brightness, and contrast.

## Getting Started

### Prerequisites

*   Julia (developed on v1.9+)

### Running the Simulations

1.  **Clone the repository.**
2.  **Navigate to the project directory.**
3.  **Activate the project environment and install dependencies:**
    ```bash
    julia --project=. -e "using Pkg; Pkg.instantiate()"
    ```

#### To Generate the Animation:

The main script for the final render is located at `scripts/8-full-render/Render.jl`.

1.  **First, run the physics computation:**
    *   Open `scripts/8-full-render/Render.jl`.
    *   Set `RUN_PHYSICS = true` and `RUN_RENDER = false`.
    *   Execute the script from the Julia REPL or command line:
        ```bash
        julia scripts/8-full-render/Render.jl
        ```
    *   This will take a significant amount of time and will create a `black_hole_data_1080p.jld2` file.

2.  **Then, render the animation:**
    *   In the same script, set `RUN_PHYSICS = false` and `RUN_RENDER = true`.
    *   You can now tweak parameters in the `render_animation` function (e.g., `EXPOSURE`, colors, intensity).
    *   Execute the script again:
        ```bash
        julia scripts/8-full-render/Render.jl
        ```
    *   This will use the data file to quickly generate the final `black_hole_1080p.gif` in the same directory.
