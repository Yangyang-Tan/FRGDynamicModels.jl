# FRGDynamicModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Yangyang-Tan.github.io/FRGDynamicModels.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Yangyang-Tan.github.io/FRGDynamicModels.jl/dev/)
[![Build Status](https://github.com/Yangyang-Tan/FRGDynamicModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Yangyang-Tan/FRGDynamicModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Yangyang-Tan/FRGDynamicModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Yangyang-Tan/FRGDynamicModels.jl)

This package provides a Julia implementation for solving the functional renormalization group (FRG) flow equations for the global effective potentials of $O(N)$ theories.

## Installation
You need to install [Julia](https://julialang.org/) first. Then you can install this package typing `] add https://github.com/Yangyang-Tan/FRGDynamicModels.jl` in the Julia REPL.

## Usage
After installing the package, you can use it by typing `using FRGDynamicModels` in the Julia REPL. 

## Examples: Solving 3d $O(4)$ model
```julia
using FRGDynamicModels
sol_lpap0 = solver_modelA(;Λ=430.0, kIR=1.0, n=4.0, Tem=141.633955
)
```