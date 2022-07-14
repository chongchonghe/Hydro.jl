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


# Wave speed estimater: min(uL - aL, uR - aR) and max(uL + aL, uR + aR)
function wave_speed_davis2(g::Grid)
    lam_minus = min.(g.vxL .- g.csL, g.vxR .- g.csR)
    lam_plus = max.(g.vxL .+ g.csL, g.vxR .+ g.csR)
    return lam_minus, lam_plus
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


""" HLL Riemann solver, 1st order. Piecewise constant construction is applied
where the left and right states on the interface are simply
    U_{i+1/2},L = U_{i}
    U_{i+1/2},R = U_{i+1}
The wave-speeds are estimated via:
    S_L = min{U_L - a_L, U_R - a_R}
    S_R = max{U_L + a_L, U_R + a_R}
"""
function hll1st(g::Grid)
    # calculate the eigenvalues
    sl = similar(g.vx)
    sr = similar(g.vx)
    for i = g.jlo-1:g.jhi
        sl[i], sr[i] = util_speed_davis2(
            g.vx[i], g.vx[i+1], g.cs[i], g.cs[i+1])
    end
    # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    alpha_plus = max.(0.0, sr)
    alpha_minus = max.(0.0, -1.0 .* sl)
    # calculate the difference of flux
    fu = calc_flux(g)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (alpha_plus[i] * fu[i, k] + alpha_minus[i] * fu[i+1, k] -
            alpha_plus[i] * alpha_minus[i] * (g.u[i+1, k] - g.u[i, k])) /
            (alpha_plus[i] + alpha_minus[i])
    end
    # calculate L(u)
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    return g.lu
end


""" Second-order HLL. """
function hll2nd(g::Grid)
    reconstruct(g)
    # the S_R and S_L of Eq. 10.47, Toro book
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


""" HLL, 1st order. Based on Chapter 7.2.1 of Springel's class notes """
function hll1st(g::Grid2d)

    F, G = calc_flux(g)

    # vx
    lam_minus, lam_plus = wave_speed_1_x(g)
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        alpha_plus[i, j] = max(0.0, lam_plus[i, j])
        alpha_minus[i, j] = max(0.0, -lam_minus[i, j])
    end
    # lam_plus = g.vx .+ g.cs
    # lam_minus = g.vx .- g.cs
    # alpha_plus = similar(g.vx)
    # alpha_minus = similar(g.vx)
    # # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    # for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
    #     alpha_plus[i, j] = max(0.0, lam_plus[i, j], lam_plus[i+1, j])
    #     alpha_minus[i, j] = max(0.0, -lam_minus[i, j], -lam_minus[i+1, j])
    # end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = (alpha_plus[i, j] * F[i, j, k] +
            alpha_minus[i, j] * F[i+1, j, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.u[i+1, j, k] - g.u[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end
    # vy
    lam_minus, lam_plus = wave_speed_1_y(g)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        alpha_plus[i, j] = max(0.0, lam_plus[i, j])
        alpha_minus[i, j] = max(0.0, -lam_minus[i, j])
    end
    # lam_plus = g.vy .+ g.cs
    # lam_minus = g.vy .- g.cs
    # alpha_plus = similar(g.vy)
    # alpha_minus = similar(g.vy)
    # for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
    #     alpha_plus[i, j] = max(0.0, lam_plus[i, j], lam_plus[i, j+1])
    #     alpha_minus[i, j] = max(0.0, -lam_minus[i, j], -lam_minus[i, j+1])
    # end
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = (alpha_plus[i, j] * G[i, j, k] +
            alpha_minus[i, j] * G[i, j+1, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.u[i, j+1, k] - g.u[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


""" Second-order HLL. This should not change g.u """
function hll2nd(g::Grid2d)

    # x component
    interpolate_x(g)
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
    interpolate_y(g)
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


function hllc(g::Grid)
    reconstruct(g)
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


function hllc(g::Grid2d)
    # x component
    interpolate_x(g)
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
    interpolate_y(g)
    # no potential sqrt error
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) 
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
    sl = similar(g.cs)
    sr = similar(g.cs)
    ss = similar(g.cs)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        # # davis2
        # sl[i, j], sr[i, j] = util_speed_davis2(
        #     g.vxL[i, j], g.vxR[i, j], g.csL[i, j], g.csR[i, j])
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

