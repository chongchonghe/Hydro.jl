include("grid2.jl")

using .Grids

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



