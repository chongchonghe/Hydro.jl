""" Time integrator """

function calc_dt(g::Grid, dt_init=1.0)
    @debug "max(g.vx) = $(maximum(abs.(g.vx)))"
    @debug "max(g.cs) = $(maximum(abs.(g.cs)))"
    @debug "max(g.rho) = $(maximum(abs.(g.rho)))"
    @debug "max(g.pressure) = $(maximum(abs.(g.pressure)))"
    C = 0.8
    themax = 0.0
    for i = g.jlo:g.jhi
        themax = max(themax, abs(g.vx[i] + g.cs[i]), abs(g.vx[i] - g.cs[i]))
    end
    return C * g.dx / themax
end

function calc_dt_new_test(g::Grid2d, dt_init)
    C = 0.9
    dt = Inf
    for j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        dt = min(dt, g.dx / max(abs(g.vx[i, j] + g.cs[i, j]), abs(g.vx[i, j] - g.cs[i, j])))
        dt = min(dt, g.dy / max(abs(g.vy[i, j] + g.cs[i, j]), abs(g.vy[i, j] - g.cs[i, j])))
    end
    return min(C * dt, 1.05 * dt_init)
end


function calc_dt(g::Grid2d, dt_init=inf)
    C = 0.8
    themax = 0.0
    # @debug "max(g.vx) = $(maximum(abs.(g.vx)))"
    # @debug "max(g.vy) = $(maximum(abs.(g.vy)))"
    # @debug "max(g.cs) = $(maximum(abs.(g.cs)))"
    for j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        vx = g.vx[i, j]
        vy = g.vy[i, j]
        cs = g.cs[i, j]
        themax = max(themax,
                     max(abs(vx + cs), abs(vx - cs)) / g.dx + max(abs(vy + cs), abs(vy - cs)) / g.dy)
    end
    return C / themax
    # return C / maximum(max.(abs.(g.vx .+ g.cs), abs.(g.vx .- g.cs)) ./ g.dx .+ max.(abs.(g.vy .+ g.cs), abs.(g.vy .- g.cs)) ./ g.dy)
end

# General Euler, 1st order
function euler(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    lu = solver(g, reconstruct)
    @. g.u = g.u + dt * lu
    rebuild(g)
    g.t += dt
    return
end

# RK2, general
function RK2(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    uold = copy(g.u)
    lu = solver(g, reconstruct)
    k1 = dt .* lu
    @. g.u = g.u + 0.5 * k1
    rebuild(g)
    lu = solver(g, reconstruct)
    k2 = dt .* lu
    @. g.u = uold + k2
    rebuild(g)
    g.t += dt
end

# RK2, general
function RK2new(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    uold = copy(g.u)
    lu = solver(g, reconstruct)
    k1 = dt .* lu
    @. g.u = g.u + k1
    rebuild(g)
    lu = solver(g, reconstruct)
    k2 = dt .* lu
    @. g.u = uold + (k1 + k2) / 2
    rebuild(g)
    g.t += dt
end

# RK3, general
function RK3(g, dt, solver::Function, reconstruct::Function, rebuild::Function)
    uold = copy(g.u)
    mat = [1. 0. 1.; 0.75 0.25 0.25; 1/3 2/3 2/3]
    for i = 1:size(mat, 1)
        lu = solver(g, reconstruct)
        # v = @view lu[:, :, 4]
        # println("lu[:, :, 4] = ", maximum(v), " ", minimum(v))
        @. g.u = mat[i, 1] * uold + mat[i, 2] * g.u + mat[i, 3] * dt * lu
        rebuild(g)
    end
    g.t += dt
end
