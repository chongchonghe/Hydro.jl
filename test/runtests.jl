include("../src/Hydro.jl")
using .Hydro

# hydro(1, 128, 0.1, "/tmp/juro/hllc-1d-sod", init_sod;
#       solver=hllc, dtout=0.01, plotit=plot_standard_sod)
# hydro(2, 128, 0.1, "/tmp/juro/hllc-2d-sod", init_sod;
#       solver=hllc, dtout=0.01, plotit=plot_standard_sod)

# hydro(1, 128, 0.1, "/tmp/juro/hll-order1-1d-sod.v2", init_sod;
#       solver=hll, order=1, dtout=0.01, plotit=plot_standard_sod)
# hydro(1, 128, 0.1, "/tmp/juro/hll-1d-sod.v2", init_sod;
#       solver=hll, order=2, dtout=0.01, plotit=plot_standard_sod)
# hydro(2, 128, 0.1, "/tmp/juro/hll-order1-2d-sod.v2", init_sod;
#       solver=hll, order=1, dtout=0.01, plotit=plot_standard_sod)
# hydro(2, 128, 0.1, "/tmp/juro/hll-2d-sod.v2", init_sod;
#       solver=hll, order=2, dtout=0.01, plotit=plot_standard_sod)

# hydro(1, 128, 0.1, "/tmp/juro/hllc-order1-1d-sod.v2", init_sod;
#       solver=hllc, order=1, dtout=0.01, plotit=plot_standard_sod)
# hydro(1, 128, 0.1, "/tmp/juro/hllc-1d-sod.v2", init_sod;
#       solver=hllc, order=2, dtout=0.01, plotit=plot_standard_sod)
# hydro(2, 128, 0.1, "/tmp/juro/hllc-order1-2d-sod.v2", init_sod;
#       solver=hllc, order=1, dtout=0.01, plotit=plot_standard_sod)
# hydro(2, 128, 0.1, "/tmp/juro/hllc-2d-sod.v2", init_sod;
#       solver=hllc, order=2, dtout=0.01, plotit=plot_standard_sod)

# hydro(1, 1024, 0.2, "/tmp/juro/hll-2d-sod-1024", init_sod;
#       solver=hll, order=2, dtout=0.01, plotit=plot_standard_sod)
# hydro(1, 1024, 0.2, "/tmp/juro/hllc-2d-sod-1024", init_sod;
#       solver=hllc, order=2, dtout=0.01, plotit=plot_standard_sod)

hydro(2, 256, 2.0, "/tmp/juro/hll-2d-KH-256", init_KH;
             solver=hll, order=2, dtout=0.005, plotit=plot_heat)

# using Hydro
# using Test

# @testset "Hydro.jl" begin
#     # Write your tests here.
#     @test hydro(1, 128, 0.1, "/tmp/juro/1d-sod", init_sod; dtout=0.01, plotit=plot_standard_sod) == 0
# end
