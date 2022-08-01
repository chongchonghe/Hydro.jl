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
    prims = zeros(xlen, 3)     # rho, vx, pressure
    rho = @view prims[:, 1]
    vx = @view prims[:, 2]
    pressure = @view prims[:, 3]
    p = pressure
    cs = zeros(xlen)
    E = similar(cs)
    epsilon = similar(cs)
    uL = similar(u) # left state
    uR = similar(u) # right state
    # fuL = similar(u) # left state
    # fuR = similar(u) # right state
    primsL = similar(prims)
    rhoL = @view primsL[:, 1]
    vxL = @view primsL[:, 2]
    pressureL = @view primsL[:, 3]
    pL = pressureL
    primsR = similar(prims)
    rhoR = @view primsR[:, 1]
    vxR = @view primsR[:, 2]
    pressureR = @view primsR[:, 3]
    pR = pressureR
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
    prims = zeros(xlen, ylen, 4)     # rho, vx, pressure
    rho = @view prims[:, :, 1]
    vx = @view prims[:, :, 2]
    vy = @view prims[:, :, 3]
    pressure = @view prims[:, :, 4]
    p = pressure
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
    pL = pressureL
    primsR = similar(prims)
    rhoR = @view primsR[:, :, 1]
    vxR = @view primsR[:, :, 2]
    vyR = @view primsR[:, :, 3]
    pressureR = @view primsR[:, :, 4]
    pR = pressureR
    csL = similar(cs)
    csR = similar(cs)
end

