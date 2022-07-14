using Hydro
using Documenter

DocMeta.setdocmeta!(Hydro, :DocTestSetup, :(using Hydro); recursive=true)

makedocs(;
    modules=[Hydro],
    authors="ChongChong He <che1234@umd.edu> and contributors",
    repo="https://github.com/chongchonghe/Hydro.jl/blob/{commit}{path}#{line}",
    sitename="Hydro.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://chongchonghe.github.io/Hydro.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/chongchonghe/Hydro.jl",
    devbranch="main",
)
