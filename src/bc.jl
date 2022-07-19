""" Fill 1D boundaries: transimissive """
function fill_trans_bc(g::Grid)
    for k = 1:3, i = 1:g.ng
        g.u[i, k] = g.u[g.jlo, k]
        g.u[g.jhi+i, k] = g.u[g.jhi, k]
    end
end


""" Fill 2D boundaries: transimissive """
function fill_trans_bc(g::Grid2d)
    for i = 1:g.ng
        g.u[i, :, :] .= g.u[g.xjlo, :, :]
        g.u[g.xjhi + i, :, :] .= g.u[g.xjhi, :, :]
        g.u[:, i, :] .= g.u[:, g.yjlo, :]
        g.u[:, g.yjhi + i, :] .= g.u[:, g.yjhi, :]
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
