using Documenter
using GalacticPotentials


makedocs(
    sitename="GalacticPotentials",
    format=Documenter.HTML(),
    modules=[GalacticPotentials],
    pages=[
        "Overview" => [
            "Getting Started" => "index.md",
            "Gravitational Potentials" => "potentials.md",
            "Example Usage" => "examples.md",
        ],
        "Reference" => "reference.md"
    ]
)

deploydocs(
    target="build",
    repo="github.com/cadojo/GalacticPotentials.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable" => "v^", "manual", "v#.#", "v#.#.#"],
)
