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


# Wave speed estimater: Roe 1981 JCP 43:357–372.
# need L and R states of rho, v, E, p
function util_speed_roe(rhoL, vL, EL, pL, rhoR, vR, ER, pR)
    uhat = (sqrt(rhoL) * vL + sqrt(rhoR) * vR) / (sqrt(rhoL) + sqrt(rhoR))
    HL = (EL + pL) / rhoL
    HR = (ER + pR) / rhoR
    Hhat = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR))
    ahat = sqrt((gamma - 1) * (Hhat - 0.5 * uhat^2))
    return uhat - ahat, uhat + ahat
end


# Wave speed estimater: min(uL - aL, uR - aR) and max(uL + aL, uR + aR)
function wave_speed_davis2(g::Grid)
    lam_minus = min.(g.velL .- g.csL, g.velR .- g.csR)
    lam_plus = max.(g.velL .+ g.csL, g.velR .+ g.csR)
    return lam_minus, lam_plus
end


""" HLL Riemann solver, 1st order. Piecewise constant construction is applied
where the left and right states on the interface are simply
    U_{i+1/2},L = u_{i}
    U_{i+1/2},R = u_{i+1}
The wave-speeds are estimated via the simple expression:
    S_L = min{U_L - a_L, U_R - a_R}
    S_R = max{U_L + a_L, U_R + a_R}
"""
function hll1st(g::Grid)
    # # calculate the eigenvalues
    # lam_minus, lam_plus = wave_speed_davis2(g)

    # calculate the eigenvalues
    lam_plus = similar(g.vel)
    lam_minus = similar(g.vel)
    for i = g.jlo-1:g.jhi
        lam_minus[i], lam_plus[i] = util_speed_davis2(
            g.vel[i], g.vel[i+1], g.cs[i], g.cs[i+1])
    end

    # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    alpha_minus = max.(0.0, -lam_minus)
    alpha_plus = max.(0.0, lam_plus)
    # calculate the difference of flux
    fu = calc_flux(g)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (alpha_plus[i] * fu[i, k] + alpha_minus[i] * fu[i+1, k]
                        - alpha_plus[i] * alpha_minus[i] *
                            (g.u[i+1, k] - g.u[i, k])) /
                            (alpha_plus[i] + alpha_minus[i])
    end
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    i = g.xmid
    @debug "\ng.t = $(g.t)"
    @debug "all the terms in fhll:"
    @debug alpha_plus[i], fu[i, 2], alpha_minus[i], fu[i+1, 2], \
        g.u[i+1, 2], g.u[i, 2]
    @debug "fhll[$(i), $(2)]: $(g.fhll[i, 2])"
    return g.lu
end


""" Second-order HLL. This should not change g.u """
function hll2nd(g::Grid)
    @assert g.ng ≥ 2

    reconstruct(g)

    # alpha_plus = max.(0.0, g.velL .+ g.csL, g.velR .+ g.csR)
    # alpha_minus = max.(0.0, -g.velL .+ g.csL, -g.velR .+ g.csR)

    lam_minus, lam_plus = wave_speed_davis2(g) 
    # # = S_R and S_L of Eq. 10.47 toro
    # lam_minus = similar(g.cs)
    # lam_plus = similar(g.cs)
    # for i = g.jlo-1:g.jhi
    #     lam_minus[i], lam_plus[i] = util_speed_davis2(
    #         g.velL[i], g.velR[i], g.csL[i], g.csR[i])
    # end

    alpha_plus = max.(0.0, lam_plus) 
    alpha_minus = max.(0.0, -1.0 .* lam_minus) 

    # save some memory by calculating fhll step by step
    calc_flux!(g.rhoL, g.velL, g.pressureL, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = alpha_plus[i] * g.fu[i, k]
    end
    calc_flux!(g.rhoR, g.velR, g.pressureR, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (g.fhll[i, k] + alpha_minus[i] * g.fu[i, k] -
            alpha_plus[i] * alpha_minus[i] * (g.uR[i, k] - g.uL[i, k])) /
            (alpha_plus[i] + alpha_minus[i])
    end
    # calculate L(u)
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end

    # debug
    @debug "\ng.t = $(g.t)"
    @debug "all the terms in fhll:"
    i = g.xmid
    k = 2
    @debug alpha_plus[i], g.fuL[i, k], alpha_minus[i], g.fuR[i, k], g.uR[i, k], g.uL[i, k]
    @debug "fhll[$(i), $(k)]: $(g.fhll[i, k])"

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

    @assert g.ng ≥ 2

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
            g.rhoR[i, j], g.vxR[i, j], g.uR[i, j, 4], g.pressureR[i, j])
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
            g.rhoR[i, j], g.vyR[i, j], g.uR[i, j, 4], g.pressureR[i, j])
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
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dx
    end

    return g.lu

end
