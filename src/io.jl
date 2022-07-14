""" References:
JLD.jl
https://stackoverflow.com/questions/41708418/iterate-through-fields-of-a-composite-type-in-julia
"""

using JLD

function grid_write(g, fo::String)
    jldopen(fo, "w") do file
        write(file, "g", g)
    end
end

function from_jld_to_grid(jld, g)
    for v in fieldnames(typeof(g))
        setfield!(g, v, getfield(jld, v))
    end
end

function grid_read(g, fi::String)
    gread = jldopen(fi, "r") do file
        read(file, "g")
    end
    from_jld_to_grid(gread, g)
    return
end
