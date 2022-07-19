# Requires grid.jl, flux.jl, reconstruction.jl

""" Solver """

# uL - aL, uR - aR
function util_speed_davis1(uL, uR, csL, csR)
    return uL - csL, uR - csR
end


# min(uL - aL, uR - aR), max(uL + aL, uR + aR)
function util_speed_davis2(uL, uR, csL, csR)
    return min(uL - csL, uR - csR), max(uL + csL, uR + csR)
end


# Wave speed estimater: Roe 1981 JCP 43:357–372.
# need L and R states of rho, v, E, p
function util_speed_roe(rhoL, vL, EL, pL, rhoR, vR, ER, pR, gamma)
    uhat = (sqrt(rhoL) * vL + sqrt(rhoR) * vR) / (sqrt(rhoL) + sqrt(rhoR))
    HL = (EL + pL) / rhoL
    HR = (ER + pR) / rhoR
    Hhat = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR))
    ahat = sqrt((gamma - 1) * (Hhat - 0.5 * uhat^2))
    return uhat - ahat, uhat + ahat
end


""" Second-order HLL. """
function hll(g::Grid, reconstruct::Function=reconstruct2nd)
    reconstruct(g)
    # update uR and uL
    prim2cons!(g.rhoL, g.vxL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.vxR, g.pressureR, g.uR, g.gamma)
    # sl and sr are S_R and S_L of Eq. 10.47, Toro book
    sl = similar(g.cs)
    sr = similar(g.cs)
    for i = g.jlo-1:g.jhi
        sl[i], sr[i] = util_speed_davis2(
            g.vxL[i], g.vxR[i], g.csL[i], g.csR[i])
    end
    # this is a simplification of Eq. 10.21, Toro book
    alpha_plus = max.(0.0, sr)
    alpha_minus = max.(0.0, -1.0 .* sl)
    # save some memory by calculating fhll step by step
    calc_flux!(g.rhoL, g.vxL, g.pressureL, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = alpha_plus[i] * g.fu[i, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.pressureR, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (g.fhll[i, k] + alpha_minus[i] * g.fu[i, k] -
            alpha_plus[i] * alpha_minus[i] * (g.uR[i, k] - g.uL[i, k])) /
            (alpha_plus[i] + alpha_minus[i])
    end
    # calculate L(u)
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    return g.lu
end


""" Second-order HLL. This should not change g.u """
function hll(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    # interpolate_x(g)
    reconstruct(g, 1)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    # alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    # alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    lam_minus = similar(g.cs)
    lam_plus = similar(g.cs)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        # # davis2
        # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
    end
    alpha_minus = max.(0.0, -1.0 .* lam_minus)
    alpha_plus = max.(0.0, lam_plus)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.xfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] + alpha_minus[i, j] * g.xfu[i, j, k]
                           - alpha_plus[i, j] * alpha_minus[i, j] *
                               (g.uR[i, j, k] - g.uL[i, j, k])) /
                               (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u) contributed by the x component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end

    # y component
    # interpolate_y(g)
    reconstruct(g, 2)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    # alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    # alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    lam_minus = similar(g.cs)
    lam_plus = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        # # davis2
        # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
        #     g.vyL[i, j], g.vyR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vyL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vyR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
    end
    alpha_minus = max.(0.0, -1.0 .* lam_minus)
    alpha_plus = max.(0.0, lam_plus)

    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.yfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] +
            alpha_minus[i, j] * g.yfu[i, j, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.uR[i, j, k] - g.uL[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


""" try SL = vxL - aL for both x and y components """
function hll_vxonly(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    # interpolate_x(g)
    reconstruct(g, 1)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    # alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    # alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    lam_minus = similar(g.cs)
    lam_plus = similar(g.cs)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        # # davis2
        # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
    end
    alpha_minus = max.(0.0, -1.0 .* lam_minus)
    alpha_plus = max.(0.0, lam_plus)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.xfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] + alpha_minus[i, j] * g.xfu[i, j, k]
                           - alpha_plus[i, j] * alpha_minus[i, j] *
                               (g.uR[i, j, k] - g.uL[i, j, k])) /
                               (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u) contributed by the x component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end

    # y component
    # interpolate_y(g)
    reconstruct(g, 2)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    # alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    # alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    lam_minus = similar(g.cs)
    lam_plus = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        # # davis2
        # lam_minus[i, j], lam_plus[i, j] = util_speed_davis2(
        #     g.vyL[i, j], g.vyR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        lam_minus[i, j], lam_plus[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
    end
    alpha_minus = max.(0.0, -1.0 .* lam_minus)
    alpha_plus = max.(0.0, lam_plus)

    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.yfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] +
            alpha_minus[i, j] * g.yfu[i, j, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.uR[i, j, k] - g.uL[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


function hllc(g::Grid, reconstruct::Function=reconstruct2nd)
    reconstruct(g)
    # update uR and uL
    prim2cons!(g.rhoL, g.vxL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.vxR, g.pressureR, g.uR, g.gamma)
    # sl and sr are S_R and S_L of Eq. 10.47, Toro book
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for i = g.jlo-1:g.jhi
        # # Davis2
        # sl[i], sr[i] = util_speed_davis2(
        #     g.vxL[i], g.vxR[i], g.csL[i], g.csR[i])
        # Roe
        sl[i], sr[i] = util_speed_roe(
            g.rhoL[i], g.vxL[i], g.uL[i, 3], g.pL[i],
            g.rhoR[i], g.vxR[i], g.uR[i, 3], g.pR[i], g.gamma)
        ss[i] = (g.pR[i] - g.pL[i] + g.rhoL[i] * g.vxL[i] *
            (sl[i] - g.vxL[i]) - g.rhoR[i] * g.vxR[i] *
            (sr[i] - g.vxR[i])) / (g.rhoL[i] *
            (sl[i] - g.vxL[i]) - g.rhoR[i] * (sr[i] - g.vxR[i]))
    end
    calc_flux!(g.rhoL, g.vxL, g.pressureL, g.gamma, g.fu)
    fuL = copy(g.fu)
    calc_flux!(g.rhoR, g.vxR, g.pressureR, g.gamma, g.fu)
    fuR = g.fu
    for k = 1:3, i = g.jlo-1:g.jhi
        fsl = ss[i] * (sl[i] * g.uL[i, k] - fuL[i, k])
        fsr = ss[i] * (sr[i] * g.uR[i, k] - fuR[i, k])
        if k == 2
            fsl += sl[i] * (g.pL[i] + g.rhoL[i] * (sl[i] - g.vxL[i]) * (ss[i] - g.vxL[i]))
            fsr += sr[i] * (g.pR[i] + g.rhoL[i] * (sr[i] - g.vxR[i]) * (ss[i] - g.vxR[i]))
        elseif k == 3
            fsl += sl[i] * (g.pL[i] + g.rhoL[i] * (sl[i] - g.vxL[i]) * (ss[i] - g.vxL[i])) * ss[i]
            fsr += sr[i] * (g.pR[i] + g.rhoL[i] * (sr[i] - g.vxR[i]) * (ss[i] - g.vxR[i])) * ss[i]
        end
        fsl /= sl[i] - ss[i]
        fsr /= sr[i] - ss[i]
        if 0 <= sl[i]
            g.fhll[i, k] = fuL[i, k]
        elseif sl[i] < 0 <= ss[i]
            g.fhll[i, k] = fsl
        elseif ss[i] < 0 <= sr[i]
            g.fhll[i, k] = fsr
        else
            g.fhll[i, k] = fuR[i, k]
        end
    end
    # calculate L(u)
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    return g.lu
end


function hllc(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    # interpolate_x(g)
    reconstruct(g, 1)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        sl[i, j], sr[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vxL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * g.vxR[i, j] *
            (sr[i, j] - g.vxR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]))
    end
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.xfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.xfu
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        fsl = ss[i, j] * (sl[i, j] * g.uL[i, j, k] - fuL[i, j, k])
        fsr = ss[i, j] * (sr[i, j] * g.uR[i, j, k] - fuR[i, j, k])
        if k == 2
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j]))
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j]))
        elseif k == 4
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j])) * ss[i, j]
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j])) * ss[i, j]
        end
        fsl /= sl[i, j] - ss[i, j]
        fsr /= sl[i, j] - ss[i, j]
        if 0 <= sl[i, j]
            g.fhll[i, j, k] = fuL[i, j, k]
        elseif sl[i, j] < 0 <= ss[i, j]
            g.fhll[i, j, k] = fsl
        elseif ss[i, j] < 0 <= sr[i, j]
            g.fhll[i, j, k] = fsr
        else
            g.fhll[i, j, k] = fuR[i, j, k]
        end
    end
    # calculate L(u)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end

    # y component
    # interpolate_y(g)
    reconstruct(g, 2)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vyL[i, j], g.vyR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        sl[i, j], sr[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vyL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vyR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vyL[i, j] *
            (sl[i, j] - g.vyL[i, j]) - g.rhoR[i, j] * g.vyR[i, j] *
            (sr[i, j] - g.vyR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vyL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vyR[i, j]))
    end
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.yfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.yfu
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        fsl = ss[i, j] * (sl[i, j] * g.uL[i, j, k] - fuL[i, j, k])
        fsr = ss[i, j] * (sr[i, j] * g.uR[i, j, k] - fuR[i, j, k])
        if k == 2
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vyL[i, j]) * (ss[i, j] - g.vyL[i, j]))
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vyR[i, j]) * (ss[i, j] - g.vyR[i, j]))
        elseif k == 4
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vyL[i, j]) * (ss[i, j] - g.vyL[i, j])) * ss[i, j]
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vyR[i, j]) * (ss[i, j] - g.vyR[i, j])) * ss[i, j]
        end
        fsl /= sl[i, j] - ss[i, j]
        fsr /= sl[i, j] - ss[i, j]
        if 0 <= sl[i, j]
            g.fhll[i, j, k] = fuL[i, j, k]
        elseif sl[i, j] < 0 <= ss[i, j]
            g.fhll[i, j, k] = fsl
        elseif ss[i, j] < 0 <= sr[i, j]
            g.fhll[i, j, k] = fsr
        else
            g.fhll[i, j, k] = fuR[i, j, k]
        end
    end
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


function hllc_vxonly(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    # interpolate_x(g)
    reconstruct(g, 1)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        sl[i, j], sr[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vxL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * g.vxR[i, j] *
            (sr[i, j] - g.vxR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]))
    end
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.xfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.xfu
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        fsl = ss[i, j] * (sl[i, j] * g.uL[i, j, k] - fuL[i, j, k])
        fsr = ss[i, j] * (sr[i, j] * g.uR[i, j, k] - fuR[i, j, k])
        if k == 2
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j]))
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j]))
        elseif k == 4
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j])) * ss[i, j]
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j])) * ss[i, j]
        end
        fsl /= sl[i, j] - ss[i, j]
        fsr /= sl[i, j] - ss[i, j]
        if 0 <= sl[i, j]
            g.fhll[i, j, k] = fuL[i, j, k]
        elseif sl[i, j] < 0 <= ss[i, j]
            g.fhll[i, j, k] = fsl
        elseif ss[i, j] < 0 <= sr[i, j]
            g.fhll[i, j, k] = fsr
        else
            g.fhll[i, j, k] = fuR[i, j, k]
        end
    end
    # calculate L(u)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end

    # y component
    # interpolate_y(g)
    reconstruct(g, 2)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vyL[i, j], g.vyR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        sl[i, j], sr[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vxL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * g.vxR[i, j] *
            (sr[i, j] - g.vxR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]))
    end
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.yfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.yfu
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        fsl = ss[i, j] * (sl[i, j] * g.uL[i, j, k] - fuL[i, j, k])
        fsr = ss[i, j] * (sr[i, j] * g.uR[i, j, k] - fuR[i, j, k])
        if k == 2
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j]))
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j]))
        elseif k == 4
            fsl += sl[i, j] * (g.pL[i, j] + g.rhoL[i, j] *
                (sl[i, j] - g.vxL[i, j]) * (ss[i, j] - g.vxL[i, j])) * ss[i, j]
            fsr += sr[i, j] * (g.pR[i, j] + g.rhoL[i, j] *
                (sr[i, j] - g.vxR[i, j]) * (ss[i, j] - g.vxR[i, j])) * ss[i, j]
        end
        fsl /= sl[i, j] - ss[i, j]
        fsr /= sl[i, j] - ss[i, j]
        if 0 <= sl[i, j]
            g.fhll[i, j, k] = fuL[i, j, k]
        elseif sl[i, j] < 0 <= ss[i, j]
            g.fhll[i, j, k] = fsl
        elseif ss[i, j] < 0 <= sr[i, j]
            g.fhll[i, j, k] = fsr
        else
            g.fhll[i, j, k] = fuR[i, j, k]
        end
    end
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


# LAX scheme, 1D
function lax(g::Grid)
    fu = calc_flux(g)
    for k = 1:3, i = g.jlo:g.jhi
        # calculate L(u)
        g.lu[i, k] = -0.5 / g.dx * (fu[i + 1, k] - fu[i - 1, k])
        # replace g.u with (g.u[i-1, k] + g.u[i+1, k]) / 2
        g.u[i, k] = 0.5 * (g.u[i-1, k] + g.u[i+1, k])
    end
    # @. g.u[g.jlo:g.jhi, :] = 0.5 * (g.u[g.jlo-1:g.jhi-1, :] + g.u[g.jlo+1:g.jhi+1, :])
    return g.lu
end


# LAX scheme, 2D
function lax(g::Grid2d)
    F, G = calc_flux(g)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        # calculate L(u)
        g.lu[i, j, k] = -0.5 / g.dx * (F[i + 1, j, k] - F[i - 1, j, k]) -
            0.5 / g.dy * (G[i, j + 1, k] - G[i, j - 1, k])
        # replace g.u with (g.u[i-1, j] + g.u[i+1, j]) / 2
        g.u[i, j, k] = 0.25 * (g.u[i-1, j, k] + g.u[i+1, j, k] +
            g.u[i, j-1, k] + g.u[i, j+1, k])
    end
    return g.lu
end


function hllc2(g::Grid2d, reconstruct::Function=reconstruct2nd)
    # x component
    reconstruct(g, 1)           # reconstruct ρ, vx, vy, p, cs
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) # update uL and uR
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.xfu)

    v = @view fuL[:, :, 4]
    @debug "fuL = $(maximum(v)), $(minimum(v))"

    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.xfu

    v = @view fuR[:, :, 4]
    @debug "fuR = $(maximum(v)), $(minimum(v))"

    # Step I: pressure estimate
    p_pvrs = @. 0.5 * (g.pL + g.pR) - 0.5 * (g.vxR - g.vxL) *
        0.5 * (g.rhoL + g.rhoR) * 0.5 * (g.csL + g.csR)
    ps = max.(0.0, p_pvrs)
    # Step II: wave speed estimates
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
        # roe
        sl[i, j], sr[i, j] = util_speed_roe(
            g.rhoL[i, j], g.vxL[i, j], g.uL[i, j, 4], g.pressureL[i, j],
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j],
            g.gamma)
        # star region
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vxL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * g.vxR[i, j] *
            (sr[i, j] - g.vxR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]))
    end

    v = @view sl[g.xjlo:g.xjhi, g.yjlo:g.yjhi]
    @debug "sl = $(maximum(v)), $(minimum(v))"
    v = @view sr[g.xjlo:g.xjhi, g.yjlo:g.yjhi]
    @debug "sr = $(maximum(v)), $(minimum(v))"
    v = @view ss[g.xjlo:g.xjhi, g.yjlo:g.yjhi]
    @debug "ss = $(maximum(v)), $(minimum(v))"
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        @assert !(ss[i, j] ≈ sl[i, j])
        @assert !(ss[i, j] ≈ sr[i, j])
        @assert !(ss[i, j] ≈ g.vxL[i, j]) "i = $(i), j = $(j), ss = $(ss[i, j]), g.vxL = $(g.vxL[i, j])"
        @assert !(ss[i, j] ≈ g.vxR[i, j]) "i = $(i), j = $(j)"
    end

    v = @view g.fhll[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 4]
    @debug "g.fhll = $(maximum(v)), $(minimum(v))"

    # Step III: HLLC ﬂux.
    uslk = zeros(4)
    usrk = zeros(4)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        usl = g.rhoL[i, j] * (sl[i, j] - g.vxL[i, j]) / (sl[i, j] - ss[i, j])
        uslk[1] = usl
        uslk[2] = usl * ss[i, j]
        uslk[3] = usl * g.vxL[i, j]
        uslk[4] = g.uL[i, j, 4] / g.rhoL[i, j] + (ss[i, j] - g.vxL[i, j]) * (ss[i, j] + g.pL[i, j] / (g.rhoL[i, j] * (ss[i, j] - g.vxL[i, j])))

        # fsl = fuL[i, j, k] + sl[i, j] * (usl - g.uL[i, j])
    
        usr = g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]) / (sr[i, j] - ss[i, j])
        usrk[1] = usr
        usrk[2] = usr * ss[i, j]
        usrk[3] = usr * g.vxR[i, j]
        usrk[4] = usr * g.uR[i, j, 4] / g.rhoR[i, j] + (ss[i, j] - g.vxR[i, j]) * (ss[i, j] + g.pR[i, j] / (g.rhoR[i, j] * (ss[i, j] - g.vxR[i, j])))

        # fsr = fuR[i, j, k] + sr[i, j] * (usr - g.uR[i, j])

        @assert !isnan(maximum(uslk))
        @assert !isnan(maximum(usrk))
        @assert !isnan(maximum(g.uL[:, :, 4]))
        @assert !isnan(maximum(g.uR[:, :, 4]))

        if 0 <= sl[i, j]
            g.fhll[i, j, :] .= fuL[i, j, :]
        elseif 0 <= ss[i, j]
            g.fhll[i, j, :] .= fuL[i, j, :] .+ sl[i, j] .* (uslk[:] .- g.uL[i, j, :])
        elseif 0 <= sr[i, j]
            g.fhll[i, j, :] .= fuR[i, j, :] .+ sr[i, j] .* (usrk[:] .- g.uR[i, j, :])
        else
            g.fhll[i, j, :] .= fuR[i, j, :]
        end
    end

    v = @view g.fhll[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 4]
    @debug "g.fhll = $(maximum(v)), $(minimum(v))"

    # calculate L(u) contributed by the x component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end
        

    # y component
    reconstruct(g, 2)
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma)
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    fuL = copy(g.yfu)
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    fuR = g.yfu
    # Step I: pressure estimate
    p_pvrs = @. 0.5 * (g.pL + g.pR) - 0.5 * (g.vxR - g.vxL) *
        0.5 * (g.rhoL + g.rhoR) * 0.5 * (g.csL + g.csR)
    ps = max.(0.0, p_pvrs)
    # Step II: wave speed estimates
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        if ps[i, j] <= g.pL[i, j]
            qL = 1.0
        else
            qL = sqrt(1 + (g.gamma + 1) / (2 * g.gamma) * (ps[i, j] / g.pL[i, j] - 1))
        end
        if ps[i, j] <= g.pR[i, j]
            qR = 1.0
        else
            qR = sqrt(1 + (g.gamma + 1) / (2 * g.gamma) * (ps[i, j] / g.pR[i, j] - 1))
        end
        sl[i, j] = g.vxL[i, j] - g.csL[i, j] * qL
        sr[i, j] = g.vxR[i, j] + g.csR[i, j] * qR
        ss[i, j] = (g.pR[i, j] - g.pL[i, j] + g.rhoL[i, j] * g.vxL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * g.vxR[i, j] *
            (sr[i, j] - g.vxR[i, j])) / (g.rhoL[i, j] *
            (sl[i, j] - g.vxL[i, j]) - g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]))
    end
    # Step III: HLLC ﬂux.
    uslk = zeros(4)
    usrk = zeros(4)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        usl = g.rhoL[i, j] * (sl[i, j] - g.vxL[i, j]) / (sl[i, j] - ss[i, j])
        uslk[1] = usl
        uslk[2] = usl * ss[i, j]
        uslk[3] = usl * g.vxL[i, j]
        uslk[4] = g.uL[i, j, 4] / g.rhoL[i, j] + (ss[i, j] - g.vxL[i, j]) * (ss[i, j] + g.pL[i, j] / (g.rhoL[i, j] * (ss[i, j] - g.vxL[i, j])))

        # fsl = fuL[i, j, k] + sl[i, j] * (usl - g.uL[i, j])
    
        usr = g.rhoR[i, j] * (sr[i, j] - g.vxR[i, j]) / (sr[i, j] - ss[i, j])
        usrk[1] = usr
        usrk[2] = usr * ss[i, j]
        usrk[3] = usr * g.vxR[i, j]
        usrk[4] = usr * g.uR[i, j, 4] / g.rhoR[i, j] + (ss[i, j] - g.vxR[i, j]) * (ss[i, j] + g.pR[i, j] / (g.rhoR[i, j] * (ss[i, j] - g.vxR[i, j])))

        # fsr = fuR[i, j, k] + sr[i, j] * (usr - g.uR[i, j])

        if 0 <= sl[i, j]
            g.fhll[i, j, :] .= fuL[i, j, :]
        elseif sl[i, j] < 0 <= ss[i, j]
            g.fhll[i, j, :] .= fuL[i, j, :] .+ sl[i, j] .* (uslk[:] .- g.uL[i, j, :])
        elseif ss[i, j] < 0 <= sr[i, j]
            g.fhll[i, j, :] .= fuR[i, j, :] .+ sr[i, j] .* (usrk[:] .- g.uR[i, j, :])
        elseif 0 > sr[i, j]
            g.fhll[i, j, :] .= fuR[i, j, :]
        end
    end
    @debug "g.fhll = $(g.fhll[10, 10, 4])"
    # calculate L(u) contributed by the y component of the flux
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end

    return g.lu
end


function hllc3(g::Grid2d, reconstruct::Function=reconstruct2nd)
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

    @. g.lu += - (fuR - fuL) / g.dy
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        # replace g.u with (g.u[i-1, j] + g.u[i+1, j]) / 2
        g.u[i, j, k] = 0.25 * (g.u[i-1, j, k] + g.u[i+1, j, k] +
            g.u[i, j-1, k] + g.u[i, j+1, k])
    end

    # # Step I: wave speed estimate: S_L and S_R
    # for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
    #     vl = g.vxL[i, j]
    #     vr = g.vxR[i, j]
    #     rhol = g.rhoL[i, j]
    #     rhor = g.rhoR[i, j]
    #     pl = g.pressureL[i, j]
    #     pr = g.pressureR[i, j]
    #     # # davis2
    #     # sl, sr = util_speed_davis2(vl, vr, g.csL[i, j], g.csR[i, j])
    #     # roe
    #     sl, sr = util_speed_roe(
    #         rhol, vl, g.uL[i, j, 4], pl,
    #         rhor, vr, g.uR[i, j, 4], pr,
    #         g.gamma)
    #     # star region
    #     ss = (pr - pl + rhol * vl * (sl - vl) - rhor * vr *
    #         (sr - vr)) / (rhol * (sl - vl) - rhor * (sr - vr))
    #     # Mean pressure in the Star Region
    #     PLR = 0.5 * (pl + pr + rhol * (sl - vl) * (ss - vl) + rhor *
    #         (sr - vr) * (ss - vr))
    #     # Ds = [0, 1, 0, ss]
    #     ds = [0.0, 1.0, 0.0, ss]
    #     if 0 <= sl
    #         g.fhll[i, j, :] .= fuL[i, j, :]
    #     elseif 0 <= ss
    #         g.fhll[i, j, :] .= (ss .* (sl .* g.uL[i, j, :] .- fuL[i, j, :]) .+
    #             sl * PLR .* ds) ./ (sl - ss)
    #     elseif 0 <= sr
    #         g.fhll[i, j, :] .= (ss .* (sr .* g.uR[i, j, :] .- fuR[i, j, :]) .+
    #             sr * PLR .* ds) ./ (sr - ss)
    #     else
    #         g.fhll[i, j, :] .= fuR[i, j, :]
    #     end
    # end
    # # calculate L(u) contributed by the y component of the flux
    # for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
    #     g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    # end
    # v = g.lu[g.xjlo:g.xjhi, g.yjlo:g.yjhi, 2]
    # @debug "g.lu[:, :, 2] = $(maximum(v)), $(minimum(v))"

    return g.lu
end
