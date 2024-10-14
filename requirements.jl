import Pkg
requirements = Set([
    "Dates",
    "Parameters",
    "Statistics",
    "LinearAlgebra",
    "HypothesisTests",
    "Random",
    "Distances",
    "JLD2",
    "Plots",
    "ProgressBars",
    "Printf",
    "LaTeXStrings",
    "StatsBase",
    "Distributions"
])
for a_pkg in requirements
    Pkg.add(a_pkg)
end

println("Packages Loaded. :) Proceed")