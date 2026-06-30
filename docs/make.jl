using Hamersley
using Documenter

DocMeta.setdocmeta!(Hamersley, :DocTestSetup, :(using Hamersley); recursive=true)

makedocs(;
    modules=[Hamersley],
    authors="Nicholas Anton Collins-Craft",
    sitename="Hamersley.jl",
    format=Documenter.HTML(;
        canonical="https://nickcollins-craft.github.io/Hamersley.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nickcollins-craft/Hamersley.jl",
    devbranch="main",
)
