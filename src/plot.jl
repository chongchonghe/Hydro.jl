include("tools/SodShockTubeExact.jl")

using .SodShockTube
using Printf
using Plots


function twod2oned(g::Grid2d, axis='x')
    if axis == 'x'
        g1 = Grid(nx=g.nx, ng=g.ng)
        midy = floor(Int16, g.ylen/2)
        @. begin
            g1.rho = g.rho[:, midy]
            g1.vx = g.vx[:, midy]
            g1.pressure = g.pressure[:, midy]
        end
    else
        g1 = Grid(nx=g.ny, ng=g.ng)
        midx = floor(Int16, g.xlen/2)
        @. begin
            g1.rho = g.rho[midx, :]
            g1.vx = g.vy[midx, :]
            g1.pressure = g.pressure[midx, :]
        end
    end
    g1.t = g.t
    prim2cons(g1)
    cons2prim(g1)
    return g1
end


function plotnone(g; fn="", is_save=true)
    return
end


function plot_curve(g::Grid; fn="t.png", is_save=true)
    x = g.x[g.jlo:g.jhi]
    data = zeros(g.nx, 4)
    data[:, 1] .= g.rho[g.jlo:g.jhi]
    data[:, 2] .= g.pressure[g.jlo:g.jhi]
    data[:, 3] .= g.vx[g.jlo:g.jhi]
    data[:, 4] .= g.epsilon[g.jlo:g.jhi]
    p = scatter(x, data, layout=4, ms=1, legend=false, xlabel="x",
                ylabel=["rho" "p" "vel" "epsilon"], xlim=[0, 1],
                ylim=[(0., 1.1) (-0., 1.2) (-.2, 1) (1.6, 3.)],
                dpi=300, title=@sprintf("t = %.04f", g.t))
    if is_save
        savefig(fn)
    else
        return p
    end
    return
end


""" Plot 1D curve from a 2d simulation """
function plot_curve(g::Grid2d; fn="t.png", is_save=true)
    g1 = twod2oned(g)
    plot_curve(g1; fn=fn, is_save=is_save)
end


function plot_standard_sod(g::Grid; fn="t.png", is_save=true)
    p = plot_curve(g; fn=fn, is_save=false)
    # get the exact solutions
    # Set up a shock tube problem
    problem = ShockTubeProblem(
        geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = g.t,
        γ = 1.4
    )
    x = g.x[g.jlo:g.jhi]
    positions, regions, values = solve(problem, x)
    exact_e = values.p ./ (g.gamma - 1) ./ values.ρ
    plot!(x, hcat(values.ρ, values.p, values.u, exact_e), layout=4,
          color=:blue)

    # find the relative difference
    ρ = g.rho[g.jlo:g.jhi]
    relerror = sqrt(sum((ρ .- values.ρ).^2)) / sum(values.ρ)
    println("Relative error on rho is ", relerror)
    if is_save
        savefig(fn)
    else
        return p
    end
    return relerror 
end


function plot_standard_sod(g::Grid2d; fn="t.png", is_save=true)
    g1 = twod2oned(g)
    return plot_standard_sod(g1; fn=fn, is_save=is_save)
end


function plot_standard_sod_y(g::Grid2d; fn="t.png", is_save=true)
    g1 = twod2oned(g, 'y')
    return plot_standard_sod(g1; fn=fn, is_save=is_save)
end


""" Plot heatmap of the density for a 2d simulation """
# function plot_heat(g::Grid2d; fn="heat.png", is_save=true, zmin=0., zmax=2.)
function plot_heat(g::Grid2d; fn="heat.png", is_save=true, zmin=0.9, zmax=2.1)
    # calculate rho, u, p e
    x = g.x[g.xjlo:g.xjhi]
    y = g.y[g.yjlo:g.yjhi]
    z = g.rho[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    thesize = (800, 800)
    p0 = heatmap(y, x, z, dpi=300, size=thesize, clim=(zmin, zmax),
                 xlim=(0, 1), ylim=(0, 1), showaxis=false, showticks=false,
                 c=:thermal, aspectratio=:equal, colorbar=false)
    if is_save
        savefig(fn)
    else
        return p0
    end
end


""" Plot curve for a 1d simulation and heat for a 2d simulation """
function plot_curve_or_heat(g::Grid; fn="t.png", is_save=true)
    plot_curve(g; fn=fn, is_save=is_save)
end


""" Plot curve for a 1d simulation and heat for a 2d simulation """
function plot_curve_or_heat(g::Grid2d; fn="t.png", is_save=true)
    plot_heat(g; fn=fn, is_save=is_save)
end


""" Plot heatmap of the density, 2D case """
function plot_heat_four_panels(g::Grid2d; fn="heat.png", is_save=true, zmin=0., zmax=2.)
    # calculate rho, u, p e
    x = g.x[g.xjlo:g.xjhi]
    y = g.y[g.yjlo:g.yjhi]
    z1 = g.rho[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z2 = g.vx[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z3 = g.pressure[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z4 = g.vy[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    # thesize = (1600, 1600)
    # p0 = heatmap(y, x, z1, dpi=300,
    #              layout=4, title=["rho", "vx", "pressure", "vy"],
    #              size=thesize,
    #              clim=[(zmin, zmax), (-2, 2), (0, 3), (-2, 2)],
    #              xlim=(0, 1), ylim=(0, 1), showaxis=false, c=:thermal,
    #              aspectratio=:equal)
    # new style
    p1 = heatmap(y, x, z1, title="rho", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(0, 1.1))
    p2 = heatmap(y, x, z2, title="vx", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(-1, 1))
    p3 = heatmap(y, x, z3, title="pressure", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(0, 2))
    p4 = heatmap(y, x, z4, title="vy", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(-1, 1))
    p0 = plot(p1, p2, p3, p4, size=(1800, 1800), dpi=300, layout=@layout [a b; c d])
    if is_save
        savefig(fn)
    else
        return p0
    end
end
