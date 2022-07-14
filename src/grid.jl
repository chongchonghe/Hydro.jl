using Parameters 

""" Memory use:
N * 5                            # cs, E, epsilon, csL, csR
3N * 9                           # u, fu, lu, fhll, uL, uR, prims, primsL, primsR
tot memory = 32N
"""
@with_kw mutable struct Grid
    nx::Int
    ng::Int
    t = 0.
    xmin = 0.
    xmax = 1.
    xlen = nx + 2 * ng
    xmid = floor(Int64, xlen/2)
    gamma = 1.4
    lambda = gamma - 1
    jlo = ng + 1
    jhi = nx + ng
    # construct the grids
    dx = (xmax - xmin) / nx
    # x = xmin .+ (1:(nx + 2*ng) .- ng) .* dx
    xl = (collect(0:xlen-1) .- ng) .* dx # left edges
    x = xl .+ dx/2                       # centers
    u = zeros(xlen, 3)
    fu = zeros(xlen, 3)
    lu = zeros(xlen, 3)
    fhll = zeros(xlen, 3)
    prims = zeros(xlen, 3)     # rho, vel, pressure
    rho = @view prims[:, 1]
    vel = @view prims[:, 2]
    pressure = @view prims[:, 3]
    cs = zeros(xlen)
    E = similar(cs)
    epsilon = similar(cs)
    uL = similar(u) # left state
    uR = similar(u) # right state
    # fuL = similar(u) # left state
    # fuR = similar(u) # right state
    primsL = similar(prims)
    rhoL = @view primsL[:, 1]
    velL = @view primsL[:, 2]
    pressureL = @view primsL[:, 3]
    primsR = similar(prims)
    rhoR = @view primsR[:, 1]
    velR = @view primsR[:, 2]
    pressureR = @view primsR[:, 3]
    csL = similar(cs)
    csR = similar(cs)
end


@with_kw mutable struct Grid2d
    nx::Int
    ny::Int
    ng::Int
    t = 0.
    xmin = 0.
    xmax = 1.
    ymin = 0.
    ymax = 1.
    xlen = nx + 2 * ng
    ylen = ny + 2 * ng
    xmid = floor(Int64, xlen/2)
    ymid = floor(Int64, ylen/2)
    gamma = 1.4
    lambda = gamma - 1
    xjlo = ng + 1
    xjhi = nx + ng
    yjlo = ng + 1
    yjhi = ny + ng
    # the grids
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    # x = xmin .+ (collect(1:xlen) .- ng) .* dx
    # y = ymin .+ (collect(1:ylen) .- ng) .* dy
    xl = (collect(0:xlen-1) .- ng) .* dx # left edges
    x = xl .+ dx/2                       # centers
    yl = (collect(0:ylen-1) .- ng) .* dy # left edges
    y = yl .+ dy/2                       # centers
    u = zeros(xlen, ylen, 4)
    xfu = similar(u)
    yfu = similar(u)
    lu = similar(u)
    fhll = similar(u)
    prims = zeros(xlen, ylen, 4)     # rho, vel, pressure
    rho = @view prims[:, :, 1]
    vx = @view prims[:, :, 2]
    vy = @view prims[:, :, 3]
    pressure = @view prims[:, :, 4]
    cs = zeros(xlen, ylen)
    E = similar(cs)
    epsilon = similar(cs)
    uL = similar(u) # left state
    uR = similar(u) # right state
    # fuL = similar(u) # left state
    # fuR = similar(u) # right state
    primsL = similar(prims)
    rhoL = @view primsL[:, :, 1]
    vxL = @view primsL[:, :, 2]
    vyL = @view primsL[:, :, 3]
    pressureL = @view primsL[:, :, 4]
    primsR = similar(prims)
    rhoR = @view primsR[:, :, 1]
    vxR = @view primsR[:, :, 2]
    vyR = @view primsR[:, :, 3]
    pressureR = @view primsR[:, :, 4]
    csL = similar(cs)
    csR = similar(cs)
end

@with_kw mutable struct Grid2dold
    nx::Int
    ny::Int
    ng::Int
    t = 0.
    xmin = 0.
    xmax = 1.
    ymin = 0.
    ymax = 1.
    xlen = nx + 2 * ng
    ylen = ny + 2 * ng
    xmid = floor(Int64, xlen/2)
    ymid = floor(Int64, ylen/2)
    gamma = 1.4
    lambda = gamma - 1
    xjlo = ng + 1
    xjhi = nx + ng
    yjlo = ng + 1
    yjhi = ny + ng
    # the grids
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)
    x = xmin .+ (collect(1:xlen) .- ng) .* dx
    y = ymin .+ (collect(1:ylen) .- ng) .* dy
    u = zeros(xlen, ylen, 4)
    xfu = similar(u)
    yfu = similar(u)
    lu = similar(u)
    fhll = similar(u)
    prims = zeros(xlen, ylen, 4)     # rho, vel, pressure
    rho = @view prims[:, :, 1]
    vx = @view prims[:, :, 2]
    vy = @view prims[:, :, 3]
    pressure = @view prims[:, :, 4]
    cs = zeros(xlen, ylen)
    E = similar(cs)
    epsilon = similar(cs)
    uL = similar(u) # left state
    uR = similar(u) # right state
    # fuL = similar(u) # left state
    # fuR = similar(u) # right state
    primsL = similar(prims)
    rhoL = @view primsL[:, :, 1]
    vxL = @view primsL[:, :, 2]
    vyL = @view primsL[:, :, 3]
    pressureL = @view primsL[:, :, 4]
    primsR = similar(prims)
    rhoR = @view primsR[:, :, 1]
    vxR = @view primsR[:, :, 2]
    vyR = @view primsR[:, :, 3]
    pressureR = @view primsR[:, :, 4]
    csL = similar(cs)
    csR = similar(cs)
end


""" Fill 1D boundaries: transimissive """
function fill_trans_bc(g::Grid)
    for k = 1:3, i = 1:g.ng
        g.u[i, k] = g.u[g.jlo, k]
        g.u[g.jhi+i, k] = g.u[g.jhi, k]
    end
end


""" Fill 2D boundaries: transimissive """
function fill_trans_bc(g::Grid2d)
    for k = 1:4, i = 1:g.ng
        g.u[i, :, k] = g.u[g.xjlo, :, k]
        g.u[g.xjhi + i, :, k] = g.u[g.xjhi, :, k]
        g.u[:, i, k] = g.u[:, g.yjlo, k]
        g.u[:, g.yjhi + i, k] = g.u[:, g.yjhi, k]
    end
end


""" Fill periodic boundary conditions in 2D space """
function fill_periodic_bc(g::Grid2d)
    for k = 1:4, i = 1:g.ng
        g.u[i, :, k] = g.u[g.xjhi - g.ng + i, :, k]
        g.u[g.xjhi + i, :, k] = g.u[g.xjlo + i - 1, :, k]
        g.u[:, i, k] = g.u[:, g.yjhi - g.ng + i, k]
        g.u[:, g.yjhi + i, k] = g.u[:, g.yjlo + i - 1, k]
    end
end


function twod2oned(g::Grid2d)
    g1 = Grid(nx=g.nx, ng=g.ng)
    g1.t = g.t
    midy = floor(Int16, g.ylen/2)
    @. begin
        g1.rho = g.rho[:, midy]
        g1.vel = g.vx[:, midy]
        g1.pressure = g.pressure[:, midy]
    end
    prim2cons(g1)
    cons2prim(g1)
    return g1
end

