using AtomsUtilityCalculators
using Documenter

DocMeta.setdocmeta!(AtomsUtilityCalculators, :DocTestSetup, :(using AtomsUtilityCalculators); recursive=true)

makedocs(;
    modules=[AtomsUtilityCalculators],
    authors="Teemu Järvinen <teemu.j.jarvinen@gmail.com> and contributors",
    sitename="AtomsUtilityCalculators.jl",
    format=Documenter.HTML(;
        canonical="https://tjjarvinen.github.io/AtomsUtilityCalculators.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tjjarvinen/AtomsUtilityCalculators.jl",
    devbranch="main",
)
