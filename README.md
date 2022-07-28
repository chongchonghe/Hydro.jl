# Hydro

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chongchonghe.github.io/Hydro.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chongchonghe.github.io/Hydro.jl/dev/)
[![Build Status](https://github.com/chongchonghe/Hydro.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/chongchonghe/Hydro.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/chongchonghe/Hydro.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/chongchonghe/Hydro.jl)

Hydro.jl is a modular hydrodynamic code written in pure Julia. 

## Tests and Demo

Click on the gif to see a full video on vimeo.com. Click [here](https://github.com/chongchonghe/shared-files/blob/main/github/hll-2d-KH-512.new3.30hz.mp4?raw=true) to download a lossless video (34 MB).

[![gif](https://videoapi-muybridge.vimeocdn.com/animated-thumbnails/image/da398ffd-45f7-42cd-ba67-034e95081831.gif?ClientID=vimeo-core-prod&Date=1659049549&Signature=9a38f58dc1c3fcad84bd2b2165e4b7b5321ce23f)](https://vimeo.com/734536881)

## Technical Description

- Riemann solver
	  - HLL
	  - HLLC
	  - Lax
- Integrator
	  - Euler
	  - RK2
	  - RK3
- Initial conditions
	  - 2D Kelvin-Helmholtz instability
	  - 1D and 2D sod shock tube
- Boundary condition
	  - Transimissive
	  - Periodic

TODO: add descriptions to the Riemann solvers, interpolation method, etc.

## Getting Started

### Dependencies

- [Julia](https://julialang.org/), a high-level, high-performance, dynamic programming language

### Installation

This package can be installed using the Julia package manager from the Julia command line:

``` julia
julia> using Pkg; Pkg.add("ArgParse")
```

<!-- First, clone this repository to a directory -->

<!--     git clone https://github.com/chongchonghe/Juro.git -->

### Usage

Then, you can use this module either with command line interface or Julia REPL.

### Run with the command line interface

Run `julia -h` for a detailed instruction. Here are some simple examples:

    julia Juro.jl/src/run.jl sod 128 0.1 tmp
    julia Juro.jl/src/run.jl KH 512 1.0 examples/KH_512

Once the run is done, you can restart from the last snapshot by adding `--restart 10` and keeping every other parameters unchanged except `tend`. A snapshot to resume the simulation from is store in `output_dir/data`. 

Optionally, you may add a `-i` option to the `julia` command to make the run more flexible by bringing two extra features:

1.  You may interrupt the run at any time by pressing `Ctrl-c` and the program will save all the necessary data needed to resume the run before quitting.
2.  After the program is finished, you will be returned to the Julia REPL where the program was running, and you will be able to start another run immediately by using `hydro(...)`. See the next section `Run with REPL` for details.

### Run with REPL

You can use Hydro.jl via Julia REPL (read-eval-print loop).

    julia> include("Hydro.jl/src/Hydro.jl")
    
    julia> using .Hydro
    
    julia> hydro(1, 512, 0.1, "tmp", init_sod)
    
    julia> hydro(1, 512, 0.2, "tmp", init_sod, restart=10)



## Help

## Authors

ChongChjong He (che1234 @ umd.edu) (https://www.astro.umd.edu/~chongchong/)

## Version History

