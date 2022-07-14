include("SodShockTubeExact.jl")

using .SodShockTube

# Set up a shock tube problem
problem = ShockTubeProblem(
    geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
    left_state = (ρ = 1.0, u = 0.0, p = 1.0),
    right_state = (ρ = 0.125, u = 0.0, p = 0.1),
    t = 0.1,
    γ = 1.4
)

# xs = LinRange(0.0, 1.0, 1024); # x locations at which to solve
# grids: mid points of 256 grids between 0 and 1
n = 2^8
dx = 1.0 / n
xs = collect(0:n-1) .* dx .+ dx/2
# println(xs)

positions, regions, values = solve(problem, xs)

println(positions)
println(regions)
println(values)
# values.x, values.ρ; # density
# values.x, values.u; # velocity
# values.x, values.p; # pressure
# values.x, values.e; # energy (= epsilon * rho)
