""" From primilary variables (rho, vel, pressure) to conservative variables """
function prim2cons!(ρ, vel, pressure, u::Matrix{Float64}, gamma::Float64)
    u[:, 1] .= ρ
    u[:, 2] .= ρ .* vel
    @. u[:, 3] = pressure / (gamma - 1) + 0.5 * ρ * vel^2
    return
end

function prim2cons(g::Grid)
    prim2cons!(g.rho, g.vel, g.pressure, g.u, g.gamma)
    return
end


""" 2d: From primilary variables (rho, vel, pressure) to conservative variables """
function prim2cons!(ρ, vx, vy, pressure, u::Array{Float64, 3}, gamma::Float64)
    u[:, :, 1] .= ρ
    u[:, :, 2] .= ρ .* vx
    u[:, :, 3] .= ρ .* vy
    @. u[:, :, 4] = pressure / (gamma - 1) + 0.5 * ρ * (vx^2 + vy^2)
    return
end

function prim2cons(g::Grid2d)
    prim2cons!(g.rho, g.vx, g.vy, g.pressure, g.u, g.gamma)
    return
end


""" From prims to sound speed """
function prim2cs(ρ::Float64, pressure::Float64, γ::Float64)
    return sqrt(γ * pressure / ρ)
end


""" Reconstruct rho, vel, E, epsilon, and pressure from g.u """
function cons2prim(g::Grid)
    g.rho .= g.u[:, 1]
    g.vel .= g.u[:, 2] ./ g.rho
    g.E .= g.u[:, 3]
    @. g.epsilon = g.E / g.rho - 0.5 * g.vel^2
    @. g.pressure = g.lambda * g.rho * g.epsilon
    @. g.cs = sqrt(g.gamma * g.pressure / g.rho)
    return nothing
end


""" 2d: Reconstruct rho, vel, E, epsilon, and pressure from g.u """
function cons2prim(g::Grid2d)
    g.rho .= g.u[:, :, 1]
    g.vx .= g.u[:, :, 2] ./ g.rho
    g.vy .= g.u[:, :, 3] ./ g.rho
    g.E .= g.u[:, :, 4]
    @. g.epsilon = g.E / g.rho - 0.5 * (g.vx^2 + g.vy^2)
    @. g.pressure = g.lambda * g.rho * g.epsilon
    @. g.cs = sqrt(g.gamma * g.pressure / g.rho)
    return nothing
end


""" Calculate flux from prims (rho, vel, pressure) """
function calc_flux!(rho, vel, pressure, gamma::Float64, fu)
    @. fu[:, 1] = rho * vel
    @. fu[:, 2] = rho * vel * vel + pressure
    @. fu[:, 3] = (pressure / (gamma - 1) + 0.5 * rho * vel^2 + pressure) * vel
end

function calc_flux(rho, vel, pressure, gamma::Float64)
    fu = Array{Float64}(undef, size(rho, 1), 3)
    calc_flux!(rho, vel, pressure, gamma, fu)
    return fu
end

function calc_flux(g::Grid)
    calc_flux!(g.rho, g.vel, g.pressure, g.gamma, g.fu)
    return g.fu
end


""" 2D: Calculate flux from prims (rho, vel, pressure) """
function calc_flux!(ρ, vx, vy, pressure, γ::Float64, xfu, yfu)
    @. begin
        xfu[:, :, 1] = ρ * vx
        xfu[:, :, 2] = ρ * vx^2 + pressure
        xfu[:, :, 3] = ρ * vx * vy
        xfu[:, :, 4] = (pressure / (γ - 1) + 0.5 * ρ * (vx^2 + vy^2) + pressure) * vx
        yfu[:, :, 1] = ρ * vy
        yfu[:, :, 2] = ρ * vx * vy
        yfu[:, :, 3] = ρ * vy^2 + pressure
        yfu[:, :, 4] = (pressure / (γ - 1) + 0.5 * ρ * (vx^2 + vy^2) + pressure) * vy
    end
end

function calc_flux(ρ, vx, vy, pressure, γ::Float64)
    xfu = Array{Float64}(undef, size(ρ, 1), 4)
    yfu = Array{Float64}(undef, size(ρ, 1), size(ρ, 2), 4)
    calc_flux!(ρ, vx, vy, pressure, γ, xfu, yfu)
    return xfu, yfu
end

function calc_flux(g::Grid2d)
    calc_flux!(g.rho, g.vx, g.vy, g.pressure, g.gamma, g.xfu, g.yfu)
    return g.xfu, g.yfu
end
