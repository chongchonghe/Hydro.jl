module SodShockTube

using NLsolve: nlsolve
using PartialFunctions
# using Documenter

export solve, ShockTubeProblem

sound_speed(γ, p, ρ) = √(γ * p / ρ)

function shock_tube_fn!(p1, p5, ρ1, ρ5, γ, p4)
    z = (p4[1] / p5 - 1)
    c1 = sound_speed(γ, p1, ρ1)
    c5 = sound_speed(γ, p5, ρ5)

    gm1 = γ - 1
    gp1 = γ + 1
    g2 = 2 * γ

    fact = gm1 / g2 * (c5 / c1) * z / √(1 + gp1 / g2 * z)
    fact = (1 - fact)^(g2 / gm1)
    return p1 * fact - p4[1]
end

function calculate_regions(pℓ, uℓ, ρℓ, pr, ur, ρr, γ=1.4)
    if pℓ < pr
        ρ1, p1, u1 = ρr, pr, ur
        ρ5, p5, u5 = ρℓ, pℓ, uℓ
    else
        ρ1, p1, u1 = ρℓ, pℓ, uℓ
        ρ5, p5, u5 = ρr, pr, ur
    end

    p4 = nlsolve(shock_tube_fn!$(p1, p5, ρ1, ρ5, γ), [p1]).zero[1]

    z = (p4 / p5 - 1)
    c5 = sound_speed(γ, p5, ρ5)

    gm1 = γ - 1
    gp1 = γ + 1

    gmfac1 = 0.5 * gm1 / γ
    gmfac2 = 0.5 * gp1 / γ

    fact = sqrt(1 + gmfac2 * z)
    u4 = c5 * z / (γ * fact)
    ρ4 = ρ5 * (1 + gmfac2 * z) / (1 + gmfac1 * z)

    w = c5 * fact
    p3 = p4
    u3 = u4
    ρ3 = ρ1 * (p3 / p1)^(1 / γ)
    region1 = (p1, ρ1, u1)
    region3 = (p3, ρ3, u3)
    region4 = (p4, ρ4, u4)
    region5 = (p5, ρ5, u5)
    return region1, region3, region4, region5, w
end

function calc_positions(pℓ, pr, region1, region3, w, xi, t, γ)
    p1, ρ1, u1 = region1
    p3, ρ3, u3 = region3

    c1 = sound_speed(γ, p1, ρ1)
    c3 = sound_speed(γ, p3, ρ3)

    if pℓ > pr
        xsh = xi + w * t
        xcd = xi + u3 * t
        xft = xi + (u3 - c3) * t
        xhd = xi - c1 * t
    else
        xsh = xi - w * t
        xcd = xi - u3 * t
        xft = xi - (u3 - c3) * t
        xhd = xi + c1 * t
    end

    return xhd, xft, xcd, xsh
end

function region_states(pℓ, pr, region1, region3, region4, region5)
    if pℓ > pr
        return Dict(
            "Region 1" => region1,
            "Region 2" => "RAREFACTION",
            "Region 3" => region3,
            "Region 4" => region4,
            "Region 5" => region5,
        )
    else
        return Dict(
            "Region 1" => region5,
            "Region 2" => region4,
            "Region 3" => region2,
            "Region 4" => "RAREFACTION",
            "Region 5" => region1,
        )
    end
end

function create_arrays(
    pℓ, pr, xℓ, xr, positions, state1, state3, state4, state5, x_arr, γ, t, xi
)
    xhd, xft, xcd, xsh = positions
    p1, ρ1, u1 = state1
    p3, ρ3, u3 = state3
    p4, ρ4, u4 = state4
    p5, ρ5, u5 = state5

    npts = length(x_arr)
    ρ = zeros(npts)
    p = zeros(npts)
    u = zeros(npts)
    c1 = sound_speed(γ, p1, ρ1)
    gm1 = γ - 1
    gp1 = γ + 1

    if pℓ > pr
        for (i, x) in enumerate(x_arr)
            if x < xhd
                ρ[i], p[i], u[i] = ρ1, p1, u1
            elseif x < xft
                u2 = 2 / gp1 * (c1 + (x - xi) / t)
                fact = 1 - 0.5 * gm1 * u2 / c1
                ρ2 = ρ1 * fact^(2 / gm1)
                p2 = p1 * fact^(2 * γ / gm1)
                ρ[i], p[i], u[i] = ρ2, p2, u2
            elseif x < xcd
                ρ[i], p[i], u[i] = ρ3, p3, u3
            elseif x < xsh
                ρ[i], p[i], u[i] = ρ4, p4, u4
            else
                ρ[i], p[i], u[i] = ρ5, p5, u5
            end
        end
    else
        for (i, x) in enumerate(x_arr)
            if x < xhd
                ρ[i], p[i], u[i] = ρ5, p5, u5
            elseif x < xft
                ρ[i], p[i], u[i] = ρ4, p4, u4
            elseif x < xcd
                ρ[i], p[i], u[i] = ρ3, p3, u3
            elseif x < xsh
                u2 = -2 / gp1 * (c1 + (xi - x) / t)
                fact = 1 + 0.5 * (γ - 1) * u2 / c1
                ρ2 = ρ1 * fact^(2 / gm1)
                p2 = p1 * fact^(2 * γ / gm1)
                ρ[i], p[i], u[i] = ρ2, p2, u2
            else
                ρ[i], p[i], u[i] = ρ1, p1, u1
            end
        end
    end
    return x_arr, p, ρ, u
end

function solve(left_state, right_state, geometry, t, γ, x_arr)
    pℓ, ρℓ, uℓ = left_state.p, left_state.ρ, left_state.u
    pr, ρr, ur = right_state.p, right_state.ρ, right_state.u
    xℓ, xr, xi = geometry

    if xℓ ≥ xr
        error("xℓ has to be less than xr")
    end

    if xℓ ≥ xr || xi ≤ xℓ
        error("xi has to be in between xℓ and xr")
    end

    region1, region3, region4, region5, w = calculate_regions(pℓ, uℓ, ρℓ, pr, ur, ρr, γ)
    regions = region_states(pℓ, pr, region1, region3, region4, region5)
    x_positions = calc_positions(pℓ, pr, region1, region3, w, xi, t, γ)

    pos_descriptions = (
        "Head of rarefaction", "Foot of rarefaction", "Contact Discontinuity", "Shock"
    )
    positions = Dict(desc => pos for (desc, pos) in zip(pos_descriptions, x_positions))

    x, p, ρ, u = create_arrays(
        pℓ, pr, xℓ, xr, x_positions, region1, region3, region4, region5, x_arr, γ, t, xi
    )
    energy = @. p / ((γ - 1)) + 0.5 * u^2
    values = (x=x, ρ=ρ, p=p, u=u, e=energy)
    return positions, regions, values
end

"""
    ShockTubeProblem
Contains the parameters of a shock tube problem

# Fields
`geometry::Tuple{Float64, Float64, Float64}` Contains the locations of the (left edge, right edge, initial shock location)

`left_state` Completely specified thermodynamic state of the left side of the discontinuity (NamedTuple of p, ρ, u)

`right_state` Completely specified thermodynamic state of the right side of the discontinuity (NamedTuple of p, ρ, u)

`t::Float64` The time at which the shock tube problem will be solved

`γ::Float64` The heat capacity ratio of the gas in the shock tube
"""
Base.@kwdef struct ShockTubeProblem
    geometry::Tuple{Float64,Float64,Float64}
    left_state
    right_state
    t::Float64
    γ::Float64
end

"""
    solve(s::ShockTubeProblem, x_arr)
Solve the given shock tube problem at the provided x locations.

# Returns

positions: A `Dictionary` which maps descriptive names of the regions to x coordinates

regions: A `Dictionary` which maps regions ("Region 1", "Region 2", etc) to thermodynamic states (ρ, p, u) in the shock tube solution

values: A `NamedTuple` (;x, ρ, p, e) containing the x coordinates, the density, pressure, and stagnation energy, respectively

# Example

```jldoctest;setup = :(using SodShockTube)
julia> problem = ShockTubeProblem(
    geometry = (0.0, 1.0, 0.5),
    left_state = (ρ = 1.0, u = 0.0, p = 1.0),
    right_state = (ρ = 0.125, u = 0.0, p = 0.1),
    t = 0.1,
    γ = 1.4
);

julia> xs = LinRange(0.0, 1.0, 500);

julia> positions, regions, values = solve(problem, xs);

julia> positions
Dict{String, Float64} with 4 entries:
  "Shock"                 => 0.850431
  "Foot of rarefaction"   => 0.485945
  "Head of rarefaction"   => 0.263357
  "Contact Discontinuity" => 0.685491

julia> regions
Dict{String, Any} with 5 entries:
  "Region 5" => (0.1, 0.125, 0.0)
  "Region 1" => (1.0, 1.0, 0.0)
  "Region 4" => (0.30313, 0.265574, 0.927453)
  "Region 3" => (0.30313, 0.426319, 0.927453)
  "Region 2" => "RAREFACTION"
"""
function solve(s::ShockTubeProblem, x_arr)
    return solve(s.left_state, s.right_state, s.geometry, s.t, s.γ, x_arr)
end

end
