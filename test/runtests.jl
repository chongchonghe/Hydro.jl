include("../src/Hydro.jl")
using .Hydro

# hydro(1, 128, 0.1, "/tmp/juro/hllc-1d-sod", init_sod;
#       solver=hllc, dtout=0.01, plotit=plot_standard_sod)
hydro(2, 128, 0.1, "/tmp/juro/hllc-2d-sod", init_sod;
      solver=hllc, dtout=0.01, plotit=plot_standard_sod)

# using Hydro
# using Test

# @testset "Hydro.jl" begin
#     # Write your tests here.
#     @test hydro(1, 128, 0.1, "/tmp/juro/1d-sod", init_sod; dtout=0.01, plotit=plot_standard_sod) == 0
# end
