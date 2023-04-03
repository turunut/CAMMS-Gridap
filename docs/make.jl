using CAMMS
using Documenter

DocMeta.setdocmeta!(CAMMS, :DocTestSetup, :(using CAMMS); recursive=true)

makedocs(;
    modules=[CAMMS],
    authors="Francesc Turon <fturon@cimne.upc.edu>",
    repo="https://github.com/gridap/CAMMS.jl/blob/{commit}{path}#{line}",
    sitename="CAMMS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gridap.github.io/CAMMS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gridap/CAMMS.jl",
    devbranch="main",
)
