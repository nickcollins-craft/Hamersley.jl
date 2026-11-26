#=
Build script for the Hamersley documentation.

Documenter renders the pages in `docs/src` (and the docstrings they reference), runs the
doctests with `using Hamersley` in scope, and deploys the result to GitHub Pages.
=#


using Hamersley
using Documenter

# Every doctest needs the package loaded, so we set that up once for the whole module.
DocMeta.setdocmeta!(Hamersley, :DocTestSetup, :(using Hamersley); recursive = true)

# We build the HTML site from the pages listed here, in this order.
makedocs(;
    modules = [Hamersley],
    authors = "Nicholas Anton Collins-Craft",
    sitename = "Hamersley.jl",
    format = Documenter.HTML(;
        canonical = "https://nickcollins-craft.github.io/Hamersley.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Architecture" => "architecture.md",
    ],
)

# We deploy only from the main branch.
deploydocs(;
    repo = "github.com/nickcollins-craft/Hamersley.jl",
    devbranch = "main",
)
