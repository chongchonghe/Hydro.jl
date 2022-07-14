#!/usr/bin/env julia

include("./Hydro.jl")
using .Hydro
using ArgParse


function arg_parser()
    s = ArgParseSettings()
    @add_arg_table s begin
        "problem"
        help = "The problem to solve, one of sod, sod2d, ball, KH"
        required = true
        "nx"
        help = "number of pixels in the x dimention"
        arg_type = Int64
        required = true
        "tend"
        help = "the time at which to end the simulation"
        arg_type = Float64
        required = true
        "folder"
        help = "the directory to store figures"
        required = true
        "solver"
        help = "the Riemann solver, one of lax, hll1st, hll2nd"
        default = "hll2nd"
        "integrator"
        help = "The integrator to use, one of euler, RK2, RK3"
        default = "RK3"
        "--plot"
        help = "The function to make plots with, one of curve, heat, heat4panel, or auto (automatically decide curve or heat)."
        default = "auto"
        "--dtout"
        help = "the time step of plotting figures"
        arg_type = Float64
        default = 0.01
        "--ny"
        help = "number of pixels in the x dimention; negative values means ny = nx"
        arg_type = Int64
        default = -1
        "--storealldata"
        help = "toggle always storing data at the same frequency of making figures; by default only the last snapshot is save to JLD data"
        action = :store_true
        "--restart"
        help = "the snapshot to start the simulation from; the JLD data file [folder]/data/hydro_[restart].jld is used to resume the simulation."
        arg_type = Int64
        default = -1
    end
    return parse_args(s)
end


function main1()

    args = arg_parser()

    if args["solver"] == "lax"
        solver = lax
    elseif args["solver"] == "hll1st"
        solver = hll1st
    elseif args["solver"] == "hll2nd"
        solver = hll2nd
    else
        println("Unknown solver $(args["solver"])")
        exit(1)
    end

    if args["integrator"] == "euler"
        integrator = euler
    elseif args["integrator"] == "RK2"
        integrator = RK2
    elseif args["integrator"] == "RK3"
        integrator = RK3
    else
        println("Unknown integrator $(args["integrator"])")
        exit(1)
    end

    if args["plot"] == "curve"
        plotit = plot_curve
    elseif args["plot"] == "heat"
        plotit = plot_heat
    elseif args["plot"] == "heat4panels"
        plotit = plot_heat_four_panels
    elseif args["plot"] == "auto"
        plotit = plot_curve_or_heat
    else
        println("Unknown plotting function $(args["plot"])")
        exit(1)
    end

    problem = args["problem"]
    if problem == "sod"
        dim = 1
        init = init_sod
        fillbc = fill_trans_bc
    elseif problem == "sod2d"
        dim = 2
        init = init_sod
        fillbc = fill_trans_bc
    elseif problem == "ball"
        dim = 2
        init = init_ball
        fillbc = fill_periodic_bc
    elseif problem == "KH"
        dim = 2
        init = init_KH
        fillbc = fill_periodic_bc
    elseif problem == "KH2"
        dim = 2
        init = init_KH2
        fillbc = fill_periodic_bc
    elseif problem == "KHrand"
        dim = 2
        init = init_KH_rand
        fillbc = fill_periodic_bc
    else
        println("Undefined problem $(problem)")
        exit(1)
    end

    hydro(dim, args["nx"], args["tend"], args["folder"], init;
          solver=solver, integrator=integrator, fillbc=fillbc, plotit=plotit,
          dtout=args["dtout"], ny=args["ny"],
          storealldata=args["storealldata"],
          restart=args["restart"], islog=true)

    return
end

main1()
