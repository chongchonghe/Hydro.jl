include("../src/Hydro.jl")
using .Hydro

prob = ARGS[1]
if prob == "sod"
    hydro(1, 128, 0.1, "/tmp/juro/sod-128", init_sod;
          solver=hllc, dtout=0.01, plotit=plot_standard_sod)
elseif prob == "sod2d"
    hydro(2, 128, 0.1, "/tmp/juro/sod2d-128", init_sod;
          solver=hllc, dtout=0.01, plotit=plot_standard_sod)
elseif prob == "sod2dy_hll"
    hydro(2, 128, 0.1, "/tmp/juro/sod2dy-hll-128", init_sod_y;
          solver=hll, dtout=0.01, plotit=plot_standard_sod_y)
elseif prob == "sod2dy"
    hydro(2, 128, 0.1, "/tmp/juro/sod2d-128", init_sod_y;
          solver=hllc, dtout=0.01, plotit=plot_standard_sod_y)
elseif prob == "sod2d_vxonly"
    # hydro(2, 128, 0.1, "/tmp/juro/sod2dy-hll-128", init_sod_y;
    #       solver=hll, dtout=0.01, plotit=plot_standard_sod_y)
    # hydro(2, 128, 0.1, "/tmp/juro/sod2d-hll_vxonly-128", init_sod;
    #       solver=hll_vxonly, dtout=0.01, plotit=plot_standard_sod)
    hydro(2, 128, 0.1, "/tmp/juro/sod2dy-hll_vxonly-128", init_sod_y;
          solver=hll_vxonly, dtout=0.01, plotit=plot_standard_sod_y)
    hydro(2, 128, 0.1, "/tmp/juro/sod2dy-hllc_vxonly-128", init_sod_y;
          solver=hllc_vxonly, dtout=0.01, plotit=plot_standard_sod_y)
elseif prob == "sod2d_hllc2"
    hydro(2, 32, 0.1, "/tmp/juro/sod2dy-hllc2-128", init_sod;
          solver=hllc2, dtout=0.01, plotit=plotnone)
elseif prob == "sod2d_hllc3"
    # hydro(2, 64, 0.1, "/tmp/juro/sod2d-hllc3-128", init_sod;
    #       solver=hllc3, dtout=0.01, plotit=plot_standard_sod)
    hydro(2, 64, 0.1, "/tmp/juro/sod2d_y-hllc3-128", init_sod_y;
          solver=hllc3, dtout=0.01, plotit=plot_standard_sod_y)
    # hydro(2, 64, 0.1, "/tmp/juro/sod2d-hllc3-128-heat", init_sod;
    #       solver=hllc3, dtout=0.01, plotit=plot_heat_four_panels)
elseif prob == "sod2d_y_hllc3"
    hydro(2, 64, 0.1, "local/test/sod2d_y-hllc3-128-heat", init_sod_y;
          solver=hllc3, dtout=0.01, plotit=plot_heat_four_panels)
elseif prob == "sod2d_hllc3_vy"
    hydro(2, 64, 0.1, "local/test/sod2d_y-hllc3-128.use-vy", init_sod_y;
          solver=hllc3_vy, dtout=0.01, plotit=plot_standard_sod_y)
elseif prob == "ball_hllc3"
    hydro(2, 64, 0.1, "local/test/ball-hllc3-64-heat", init_ball;
          solver=hllc3, dtout=0.01, plotit=plot_heat_four_panels)
elseif prob == "sod2dy_hllc2"
    hydro(2, 64, 0.1, "/tmp/juro/sod2dy-hllc2-128", init_sod_y;
          solver=hllc2, dtout=0.01, plotit=plot_standard_sod_y)
elseif prob == "ball"
    hydro(2, 128, 0.1, "/tmp/juro/ball-hll-128", init_ball;
          solver=hll, dtout=0.01, plotit=plot_heat)
    hydro(2, 128, 0.1, "/tmp/juro/ball-hllc-128", init_ball;
          solver=hllc, dtout=0.01, plotit=plot_heat)
elseif prob == "KH"
    hydro(2, 256, 2.0, "/tmp/juro/KH-256", init_KH;
          solver=hllc, dtout=0.005, plotit=plot_heat)
end
