using Hydro
using Test

@testset "Hydro.jl" begin
    # Write your tests here.
    @test hydro(1, 128, 0.1, "/tmp/juro/1d-sod", init_sod; dtout=0.01, plotit=plot_standard_sod) == 0
end
