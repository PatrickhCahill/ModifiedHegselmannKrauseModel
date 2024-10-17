using Parameters
using LinearAlgebra
using JLD2
using Plots
using ProgressBars
using Printf
using HypothesisTests
using Distributions
using Random
using LaTeXStrings
using StatsBase

foldername = "***_run_1" # Replace with foldername corresponding to party base data.

Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

actual_std = []
my_pred_simple = []
my_pred_adv = []
VCluster1 = nothing
for t = 1:5000

    V = Vhistory[:,t]
    P = Phistory[:,t]

    P1 = P[1]
    VCluster1= V[P1-0.15 .< V .< P1+0.15]
    push!(actual_std,std(VCluster1))


    push!(my_pred_adv,0.03 / sqrt(2) * (length(VCluster1)/1000 * 1 + 1/3 * 1)^(-0.5))
    push!(my_pred_simple,0.03 / sqrt(2) * (length(1/3 * 1 + 1/3 * 1)^(-0.5)))
end

p = plot(LinRange(0, 100, 5000), actual_std,label = "Actual " * L" \sigma", color="#3188CE", xlabel="t", size=(600, 300), margin=5Plots.mm, dpi=300, ylims=[0,maximum(actual_std)])
plot!(p, LinRange(0, 100, 5000), my_pred_adv,label = "Predicted " * L" \sigma",legend=:none)
savefig(p, "sims/$foldername/Figure 2a cluster sizes.png")

