using AtomsCalculatorsUtilities
using Documenter

DocMeta.setdocmeta!(AtomsCalculatorsUtilities, :DocTestSetup, :(using AtomsCalculatorsUtilities); recursive=true)

makedocs(;
    modules=[AtomsCalculatorsUtilities],
    authors="JuliaMolSim contributors",
    sitename="AtomsCalculatorsUtilities.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        canonical="https://JuliaMolSim.github.io/AtomsCalculatorsUtilities.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "i-PI Driver" => "ipi.md"
    ],
)

deploydocs(;
    repo="github.com/tjjarvinen/AtomsCalculatorsUtilities.jl",
    devbranch="main",
)
