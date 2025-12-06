module Constants
    export G, c, M_sun, year, AU

    # Define constants. These are all float64 as Julia default.
    const G = 6.67430e-11       # Nm^2/kg^2:    universal gravitational constant 
    const c = 2.99792458e8      # m/s:          speed of light in a vacuum 
    const M_sun = 1.989e30      # kg:           1 solar mass (mass of our sun) 
    const year = 31556952.0     # seconds:      num seconds in year
    const AU = 1.496e11         # meters:       1 astronomical unit (distance from Earth to Suns)
end