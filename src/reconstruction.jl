# Requires grid.jl, flux.jl

function minmod(x, y, z)
    0.25 * abs(sign(x) + sign(y)) * (sign(x) + sign(z)) *
        min(abs(x), abs(y), abs(z))
end


""" Interpolate prims (rhoL, rhoR, velL, velR, pressureL, pressureR),
which are used to update the cons (uL, uR) """
function reconstruct(g::Grid, theta::Float64=1.5)
    cdiff = similar(g.prims)
    # NOTE: g.rhoL, g.velL, and g.pressureL are links to the rows in
    # g.primsL, so they update along with it.
    for k = 1:3
        c = @view g.prims[:, k]     # c = rho, vel, pressure for k = 1, 2, 3
        for i = g.jlo-1:g.jhi+1
            cdiff[i, k] = 0.5 * minmod(theta * (c[i] - c[i-1]),
                                       0.5 * (c[i+1] - c[i-1]),
                                       theta * (c[i+1] - c[i]))
        end
    end
    for k = 1:3, i = g.jlo-1:g.jhi
        g.primsL[i, k] = g.prims[i, k] + cdiff[i, k]
        g.primsR[i, k] = g.prims[i+1, k] - cdiff[i+1, k]
    end
    for i = g.jlo-1:g.jhi
        g.csL[i] = prim2cs(g.rhoL[i], g.pressureL[i], g.gamma)
        g.csR[i] = prim2cs(g.rhoR[i], g.pressureR[i], g.gamma)
    end

    # update uR and uL
    prim2cons!(g.rhoL, g.velL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.velR, g.pressureR, g.uR, g.gamma)
end


# Interpolate primitive variables in the x component
function interpolate_x(g::Grid2d)
    theta = 1.5
    cdiff = similar(g.prims)
    for k = 1:4
        c = @view g.prims[:, :, k]     # c = rho, vx, vy, pressure for k = 1, 2, 3, 4
        for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi+1
            cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i-1, j]),
                                          0.5 * (c[i+1, j] - c[i-1, j]),
                                          theta * (c[i+1, j] - c[i, j]))
        end
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.primsL[i, j, k] = g.prims[i, j, k] + cdiff[i, j, k]
        g.primsR[i, j, k] = g.prims[i+1, j, k] - cdiff[i+1, j, k]
    end
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.csL[i, j] = prim2cs(g.rhoL[i, j], g.pressureL[i, j], g.gamma)
        g.csR[i, j] = prim2cs(g.rhoR[i, j], g.pressureR[i, j], g.gamma)
    end
end


function interpolate_y(g::Grid2d)
    theta = 1.5
    cdiff = similar(g.prims)
    for k = 1:4
        c = @view g.prims[:, :, k]     # c = rho, vx, vy, pressure for k = 1, 2, 3, 4
        for j = g.yjlo-1:g.yjhi+1, i = g.xjlo:g.xjhi
            cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i, j-1]),
                                          0.5 * (c[i, j+1] - c[i, j-1]),
                                          theta * (c[i, j+1] - c[i, j]))
        end
    end
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.primsL[i, j, k] = g.prims[i, j, k] + cdiff[i, j, k]
        g.primsR[i, j, k] = g.prims[i, j+1, k] - cdiff[i, j+1, k]
    end
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.csL[i, j] = prim2cs(g.rhoL[i, j], g.pressureL[i, j], g.gamma)
        g.csR[i, j] = prim2cs(g.rhoR[i, j], g.pressureR[i, j], g.gamma)
    end
end