using Hydro
using Test

@testset "Hydro.jl" begin
    # Write your tests here.
    @test hydro(1, 128, 0.1, "figures/1d-sod", init_sod; dtout=0.01, plotit=plot_standard_sod) == 0
    # @test hydro(2, 128, 0.1, "figures/2d-sod", init_sod; dtout=0.01, plotit=plot_standard_sod) == 0
    # @test hydro(2, 128, 0.1, "figures/2d-sod-y", init_sod_y; dtout=0.01, plotit=plot_standard_sod_y) == 0
end
