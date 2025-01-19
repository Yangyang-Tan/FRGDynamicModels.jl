using FRGDynamicModels
using Documenter

DocMeta.setdocmeta!(FRGDynamicModels, :DocTestSetup, :(using FRGDynamicModels); recursive=true)

makedocs(;
    modules=[FRGDynamicModels],
    authors="Yang-yang Tan",
    sitename="FRGDynamicModels.jl",
    format=Documenter.HTML(;
        canonical="https://Yangyang-Tan.github.io/FRGDynamicModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Yangyang-Tan/FRGDynamicModels.jl",
    devbranch="main",
)
