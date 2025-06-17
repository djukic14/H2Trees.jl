using H2Trees
using Documenter

DocMeta.setdocmeta!(H2Trees, :DocTestSetup, :(using H2Trees); recursive=true)

makedocs(;
    modules=[H2Trees],
    authors="djukic14 <danijel.jukic14@gmail.com> and contributors",
    sitename="H2Trees.jl",
    format=Documenter.HTML(;
        canonical="https://djukic14.github.io/H2Trees.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/djukic14/H2Trees.jl",
    devbranch="main",
)
