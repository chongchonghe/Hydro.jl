function hllc(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    reconstruct(g, 1)           # reconstruct ρ, vx, vy, p, cs

    v = g.pL[g.xjlo:g.xjhi, floor(Int, g.ny/2)]
    @debug "g.pL = $(maximum(v)), $(minimum(v))"
    v = g.pR[g.xjlo:g.xjhi, floor(Int, g.ny/2)]
    @debug "g.pR = $(maximum(v)), $(minimum(v))"

    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) # update uL and uR
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.xfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.xfu

    # @. g.lu = - (fuR - fuL) / g.dx

    # Step I: wave speed estimate: S_L and S_R
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        vxl = g.vxL[i, j]
        vxr = g.vxR[i, j]
        rhol = g.rhoL[i, j]
        rhor = g.rhoR[i, j]
        pl = g.pressureL[i, j]
        pr = g.pressureR[i, j]
        # # davis2
        # sl, sr = util_speed_davis2(vxl, vxr, g.csL[i, j], g.csR[i, j])
        # roe
        sl, sr = util_speed_roe(
            rhol, vxl, g.uL[i, j, 4], pl,
            rhor, vxr, g.uR[i, j, 4], pr,
            g.gamma)
        # star region
        ss = (pr - pl + rhol * vxl *
            (sl - vxl) - rhor * vxr *
            (sr - vxr)) / (rhol *
            (sl - vxl) - rhor * (sr - vxr))
        # Mean pressure in the Star Region
        PLR = 0.5 * (pl + pr + rhol *
            (sl - vxl) * (ss - vxl) + rhor *
            (sr - vxr) * (ss - vxr))
        # Ds = [0, 1, 0, ss]
        ds = [0.0, 1.0, 0.0, ss]
        if 0 <= sl
            g.fhll[i, j, :] .= fuL[i, j, :]
        elseif 0 <= ss
            g.fhll[i, j, :] .= (ss .* (sl .* g.uL[i, j, :] .- fuL[i, j, :]) .+
                sl * PLR .* ds) ./ (sl - ss)
        elseif 0 <= sr
            g.fhll[i, j, :] .= (ss .* (sr .* g.uR[i, j, :] .- fuR[i, j, :]) .+
                sr * PLR .* ds) ./ (sr - ss)
        else
            g.fhll[i, j, :] .= fuR[i, j, :]
        end
    end
    # calculate L(u) contributed by the x component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end
    v = g.lu[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 2]
    @debug "g.lu[:, :, 2] = $(maximum(v)), $(minimum(v))"

    # y component
    reconstruct(g, 2)           # reconstruct ρ, vx, vy, p, cs
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) # update uL and uR
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.yfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.yfu
    v = fuR[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 2]
    @debug "fuR[:, :, 2] = $(maximum(v)), $(minimum(v))"

    # @. g.lu += - (fuR - fuL) / g.dy
    # for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
    #     # replace g.u with (g.u[i-1, j] + g.u[i+1, j]) / 2
    #     g.u[i, j, k] = 0.25 * (g.u[i-1, j, k] + g.u[i+1, j, k] +
    #         g.u[i, j-1, k] + g.u[i, j+1, k])
    # end

    # Step I: wave speed estimate: S_L and S_R
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        vl = g.vxL[i, j]
        vr = g.vxR[i, j]
        rhol = g.rhoL[i, j]
        rhor = g.rhoR[i, j]
        pl = g.pressureL[i, j]
        pr = g.pressureR[i, j]
        # # davis2
        # sl, sr = util_speed_davis2(vl, vr, g.csL[i, j], g.csR[i, j])
        # roe
        sl, sr = util_speed_roe(
            rhol, vl, g.uL[i, j, 4], pl,
            rhor, vr, g.uR[i, j, 4], pr,
            g.gamma)
        # star region
        ss = (pr - pl + rhol * vl * (sl - vl) - rhor * vr *
            (sr - vr)) / (rhol * (sl - vl) - rhor * (sr - vr))
        # Mean pressure in the Star Region
        PLR = 0.5 * (pl + pr + rhol * (sl - vl) * (ss - vl) + rhor *
            (sr - vr) * (ss - vr))
        # Ds = [0, 1, 0, ss]
        ds = [0.0, 1.0, 0.0, ss]
        if 0 <= sl
            g.fhll[i, j, :] .= fuL[i, j, :]
        elseif 0 <= ss
            g.fhll[i, j, :] .= (ss .* (sl .* g.uL[i, j, :] .- fuL[i, j, :]) .+
                sl * PLR .* ds) ./ (sl - ss)
        elseif 0 <= sr
            g.fhll[i, j, :] .= (ss .* (sr .* g.uR[i, j, :] .- fuR[i, j, :]) .+
                sr * PLR .* ds) ./ (sr - ss)
        else
            g.fhll[i, j, :] .= fuR[i, j, :]
        end
    end
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    v = g.lu[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 2]
    @debug "g.lu[:, :, 2] = $(maximum(v)), $(minimum(v))"

    return g.lu
end
