using AtomsCalculatorsUtilities
using Documenter

DocMeta.setdocmeta!(AtomsCalculatorsUtilities, :DocTestSetup, :(using AtomsCalculatorsUtilities); recursive=true)

makedocs(;
    modules=[AtomsCalculatorsUtilities],
    authors="Teemu Järvinen <teemu.j.jarvinen@gmail.com> and contributors",
    sitename="AtomsCalculatorsUtilities.jl",
    format=Documenter.HTML(;
        canonical="https://tjjarvinen.github.io/AtomsCalculatorsUtilities.jl",
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
