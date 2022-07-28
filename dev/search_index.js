var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Hydro","category":"page"},{"location":"#Hydro","page":"Home","title":"Hydro","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Hydro.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Hydro]","category":"page"},{"location":"#Hydro.Grid","page":"Home","title":"Hydro.Grid","text":"Memory use: N * 5                            # cs, E, epsilon, csL, csR 3N * 9                           # u, fu, lu, fhll, uL, uR, prims, primsL, primsR tot memory = 32N\n\n\n\n\n\n","category":"type"},{"location":"#Hydro.calc_flux!-Tuple{Any, Any, Any, Any, Float64, Any, Any}","page":"Home","title":"Hydro.calc_flux!","text":"2D: Calculate flux from prims (rho, vel, pressure) \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.calc_flux!-Tuple{Any, Any, Any, Float64, Any}","page":"Home","title":"Hydro.calc_flux!","text":"Calculate flux from prims (rho, vel, pressure) \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.cons2prim-Tuple{Grid2d}","page":"Home","title":"Hydro.cons2prim","text":"2d: Reconstruct rho, vel, E, epsilon, and pressure from g.u \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.cons2prim-Tuple{Grid}","page":"Home","title":"Hydro.cons2prim","text":"Reconstruct rho, vel, E, epsilon, and pressure from g.u \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.fill_periodic_bc-Tuple{Grid2d}","page":"Home","title":"Hydro.fill_periodic_bc","text":"Fill periodic boundary conditions in 2D space \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.fill_trans_bc-Tuple{Grid2d}","page":"Home","title":"Hydro.fill_trans_bc","text":"Fill 2D boundaries: transimissive \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.fill_trans_bc-Tuple{Grid}","page":"Home","title":"Hydro.fill_trans_bc","text":"Fill 1D boundaries: transimissive \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.hll","page":"Home","title":"Hydro.hll","text":"Second-order HLL. This should not change g.u \n\n\n\n\n\n","category":"function"},{"location":"#Hydro.hll-2","page":"Home","title":"Hydro.hll","text":"Second-order HLL. \n\n\n\n\n\n","category":"function"},{"location":"#Hydro.hllc","page":"Home","title":"Hydro.hllc","text":"HLLC scheme in 1D \n\n\n\n\n\n","category":"function"},{"location":"#Hydro.hllc_var1","page":"Home","title":"Hydro.hllc_var1","text":"HLLC solver. Work in progress. Some bugs remains to be fixed. \n\n\n\n\n\n","category":"function"},{"location":"#Hydro.hydro-Tuple{Any, Any, Any, String, Function}","page":"Home","title":"Hydro.hydro","text":"A general-purpose hydro solver \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_KH-Tuple{Grid2d}","page":"Home","title":"Hydro.init_KH","text":"Kelvin-Helmholtz instability \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_KH2-Tuple{Grid2d}","page":"Home","title":"Hydro.init_KH2","text":"Kelvin-Helmholtz instability \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_ball-Tuple{Grid2d}","page":"Home","title":"Hydro.init_ball","text":"Set the 2d sod shock tube: a ball of radius 0.2 at center with density 1.0, pressure 1.0. Outside has density 0.125 and pressure 0.1 \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_sod-Tuple{Grid2d}","page":"Home","title":"Hydro.init_sod","text":"Set the 2d initial conditions: a 2D version of the 1D sod shocktube \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_sod-Tuple{Grid}","page":"Home","title":"Hydro.init_sod","text":"set the initial conditions of a sod shock tube \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.init_sod_y-Tuple{Grid2d}","page":"Home","title":"Hydro.init_sod_y","text":"Set the 2d initial conditions: a 2D version of the 1D sod shocktube in the y dimension \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.plot_curve-Tuple{Grid2d}","page":"Home","title":"Hydro.plot_curve","text":"Plot 1D curve from a 2d simulation \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.plot_curve_or_heat-Tuple{Grid2d}","page":"Home","title":"Hydro.plot_curve_or_heat","text":"Plot curve for a 1d simulation and heat for a 2d simulation \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.plot_curve_or_heat-Tuple{Grid}","page":"Home","title":"Hydro.plot_curve_or_heat","text":"Plot curve for a 1d simulation and heat for a 2d simulation \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.plot_heat_four_panels-Tuple{Grid2d}","page":"Home","title":"Hydro.plot_heat_four_panels","text":"Plot heatmap of the density, 2D case \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.prim2cons!-Tuple{Any, Any, Any, Any, Array{Float64, 3}, Float64}","page":"Home","title":"Hydro.prim2cons!","text":"2d: From primilary variables (rho, vel, pressure) to conservative variables \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.prim2cons!-Tuple{Any, Any, Any, Matrix{Float64}, Float64}","page":"Home","title":"Hydro.prim2cons!","text":"From primilary variables (rho, vel, pressure) to conservative variables \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.prim2cs-Tuple{Float64, Float64, Float64}","page":"Home","title":"Hydro.prim2cs","text":"From prims to sound speed \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.reconstruct2nd-Tuple{Grid2d, Int64}","page":"Home","title":"Hydro.reconstruct2nd","text":"Interpolate primitive variables (ρ, vx, vy, p, cs) in the x and y components \n\n\n\n\n\n","category":"method"},{"location":"#Hydro.reconstruct2nd-Tuple{Grid}","page":"Home","title":"Hydro.reconstruct2nd","text":"Interpolate prims (rhoL, rhoR, vxL, vxR, pressureL, pressureR), which are then used to update the cons (uL, uR) \n\n\n\n\n\n","category":"method"}]
}
