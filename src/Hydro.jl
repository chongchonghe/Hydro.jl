#!/usr/bin/env julia

module Hydro

export hydro, Grid, Grid2d, lax, hll, hll2nd, hllc, euler, RK2, RK3,
    init_sod, init_ball, init_KH, fill_trans_bc, fill_periodic_bc,
    plot_curve, plot_heat, plot_curve_or_heat, plot_heat_four_panels,
    plot_standard_sod


include("grid.jl")
include("io.jl")
include("flux.jl")
include("bc.jl")
include("init.jl")
include("reconstruction.jl")
include("solver.jl")
include("integrator.jl")
include("plot.jl")


""" A general-purpose hydro solver """
function hydro(dim, nx, tend, folder::String, init::Function;
               solver::Function=hll2nd, 
               order::Int=2,
               integrator::Function=RK3,
               fillbc::Function=fill_trans_bc,
               plotit::Function=plot_curve_or_heat,
               dtout::Float64=0.01, ny::Int64=-1,
               storealldata::Bool=false, restart=-1,
               islog::Bool=true)

    msg = """
Running hydro.jl, a modular 1- and 2-dimensional hydrodynamic
code written in pure Julia.
Author: Chong-Chong He (che1234@umd.edu)

Solving the following problem with parameters:
Initial condition: $(init)
nx = $(nx)
dtout = $(dtout)
tend = $(tend)
Riemann solver: $(solver)
Integrator: $(integrator)
Boundary condition: $(fillbc)
Plotting function: $(plotit)
"""
    print(msg)
    run(`mkdir -p $folder`)
    if islog
        logfile = "$(folder)/log.o"
        if isfile(logfile)
            islog = false
            for i = 1:1000
                logfile = "$(folder)/log.o$(i)"
                if !isfile(logfile)
                    islog = true
                    break
                end
            end
        end
    end
    if islog
        fw = open(logfile, "w")
        write(fw, msg)
    end

    # definition and initialization
    # if solver == lax || solver == hll1st
    #     ng = 1
    # elseif solver == hll2nd
    #     ng = 2
    # end
    ng = 2
    if dim == 1
        g = Grid(nx=nx, ng=ng)
    elseif dim == 2
        if ny < 0
            ny = nx
        end
        g = Grid2d(nx=nx, ny=ny, ng=ng)
    end
    if order == 1
        reconstruct = reconstruct1st
    elseif order == 2
        reconstruct = reconstruct2nd
    end
    function rebuild(g)
        fillbc(g)
        cons2prim(g)
        return
    end

    datadir = "$(folder)/data"
    if restart < 0
        init(g)         # fill I.C. for prims (without boundary)
        prim2cons(g)
        fillbc(g)
        cons2prim(g)
        fcount = 0
        run(`mkdir -p $(datadir)`)
        if storealldata
            fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
            grid_write(g, fn)
        end
    else
        # @assert datadir != "none"
        fcount = restart
        fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
        grid_read(g, fn)
        # println(typeof(g), g.x)
    end

    # plot the first snapshot
    error = plotit(g, fn="$(folder)/hydro_$(lpad(fcount, 5, '0')).png")
    count = 0
    msg = "\ncount = $(count), fcount = $(fcount), t = $(g.t)"
    if plotit == plot_standard_sod
        msg = msg * ", relative error = $(error)"
    end
    print(msg * "\n")
    if islog
        write(fw, msg * "\n")
    end

    # evolve and plot snapshots
    tout = g.t + dtout
    try
        while true
            dt = calc_dt(g)
            # println("t = $(g.t), dt = $(dt)")
            integrator(g, dt, solver, reconstruct, rebuild)
            @debug "count = $(count), t = $(g.t), tout = $(tout)"
            # @debug "max,min(g.vx) = $(maximum(abs.(g.vx))), $(minimum(abs.(g.vx)))"
            count += 1
            if g.t > tout || g.t ≈ tout
                fcount += 1
                # plotit(g, fn="$(folder)/f-nx$(nx)-$(lpad(fcount, 5, '0')).png")
                error = plotit(g, fn="$(folder)/hydro_$(lpad(fcount, 5, '0')).png")
                msg = "count = $(count), fcount = $(fcount), " *
                        "t = $(g.t), tout = $(tout)"
                if plotit == plot_standard_sod
                    msg = msg * ", relative error = $(error)"
                end
                print(msg * "\n")
                if islog
                    write(fw, msg * "\n")
                end
                if storealldata || (g.t > tend || g.t ≈ tend)
                    # write jld data
                    fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
                    grid_write(g, fn)
                end
                if g.t > tend || g.t ≈ tend
                    break
                end
                tout += dtout
            end
        end
    catch err
        if isa(err, InterruptException)
            fn = "$(datadir)/hydro_$(lpad(fcount, 5, '0')).jld"
            msg = "Interrupted, saving $(fn)\n"
            print(msg)
            grid_write(g, fn)
            if islog
                write(fw, msg)
            end
            # return
        else
            rethrow()
        end
    end
    println("Simulation done.")
    if islog
        write(fw, "Simulation done.\n")
        close(fw)
    end
    return 0
end

end
