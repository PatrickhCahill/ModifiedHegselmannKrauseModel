# This file will run the pipeline for a given foldername, and save the results in a new folder called sims/foldername. The pipeline is as follows:
# 1. Create the intial conditions file in the data/foldername directory.
# 2. Run the simulation and save the results in the data/foldername directory.
# 3. Run the plotting script and save the results in the sims/foldername directory.
################################ Imports ################################
using Parameters
using LinearAlgebra
using JLD2
using Plots
using ProgressBars
using Printf
using HypothesisTests
using Distributions
using Random
using Distances
using LaTeXStrings
Random.seed!(3) # Set seed

push!(LOAD_PATH, joinpath(@__DIR__, "_lib")) # Add the _lib folder to the LOAD_PATH
using model_helpers
using plotting_helpers

################################ Instantiate struct ################################

@with_kw struct Params
    dimension::Int64# The dimension of the model | Should be 1 or 2
    Nv::Int64 # Number of voters
    Np::Int64 # Number of parties
    Rvv::Float64  # The threshold for a voter consider anothers' opinion
    μvv::Float64   # The rate at which voters change their minds, called the convergence rate
    Rpv::Float64   # The threshold for a voter consider anothers' opinion
    μpv::Float64   # The rate at which voters change their minds, called the convergence rate
    σv::Float64   # The standard deviation of the noise added to the voters
    σc::Float32 = sqrt(μpv / (2 * pi^2) + μvv * (4 / 3) * Rvv^3)
    Rvp::Float64   # The threshold for a voter consider anothers' opinion
    μvp::Float64  # The rate at which voters change their minds, called the convergence rate
    Rpp::Float64  # The threshold for a voter consider anothers' opinion
    μpp::Float64   # The rate at which voters change their minds, calledx the convergence rate
    σp::Float64  # The standard deviation of the noise added to the voters

    tot_time::Float64  # Total time to run the model
    Δt::Float64 # The time step
    tot_steps::Int64 = Int64(tot_time / Δt)# The total number of time steps
    plot_step::Int64  # The period between which the plots are produced
    save_step::Int64 # How often to save the model
    ensembleInitialNoise::Float64  # The initial noise of the ensemble
    poll_noise_std::Float64 # The variance of the noise in the poll
    poll_noise_var::Float64 = (poll_noise_std * Nv)^2
    n_members::Int64
    R::Matrix{Float64} = Matrix{Float64}(I, Np + 2, Np + 2) * poll_noise_var
    inflation_factor::Float64
    poll_freq::Int64

    boundary_conditions::String # Defines which boundaryconditions
end

println("HANDLING Figure 4 and 6a alternative")
#region Figure 4 and 6a alternative
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=500, # Number of voters
    Np=4, # Number of parties
    Rvv=0.05,# The threshold for a voter consider anothers' opinion
    μvv=0.3, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.1, # The threshold for a voter consider anothers' opinion
    μpv=0.6,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0.01, # The standard deviation of the noise added to the voters
    Rvp=0.4, # The threshold for a voter consider anothers' opinion
    μvp=0.01, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.1, # The threshold for a voter consider anothers' opinion
    μpp=0.01, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.002, # The standard deviation of the noise added to the voters
    tot_time=500, # Total time to run the model
    Δt=0.05, # The time step
    plot_step=10, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="periodic"
)
foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.2 0.35 0.6 0.78], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.

(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the data from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:2:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.9, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)

orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)



if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end

savefig(p1, "sims/$(foldername)/noisyconvergence.png")
savefig(p2, "sims/$(foldername)/noisyconvergence_Dhat.png")
savefig(p3, "sims/$(foldername)/noisyconvergence_Qvv.png")



open("sims/$(foldername)/test.txt", "w") do file
    write(file, "Figure2.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 5 and 6b")
#region Figure 5 and 6b
the_params = Params(
    dimension=2, # The dimension of the model | Should be 1 or 2
    Nv=500, # Number of voters
    Np=5, # Number of parties
    Rvv=0.2,# The threshold for a voter consider anothers' opinion
    μvv=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.1, # The threshold for a voter consider anothers' opinion
    μpv=1,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0.02, # The standard deviation of the noise added to the voters
    Rvp=0.5, # The threshold for a voter consider anothers' opinion
    μvp=0.05, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.05, # The threshold for a voter consider anothers' opinion
    μpp=0.02, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.002, # The standard deviation of the noise added to the voters
    tot_time=500, # Total time to run the model
    Δt=0.1, # The time step
    plot_step=10, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="periodic")
foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 2), [0.2 0.2; 0.75 0.75; 0.6 0.8; 0.5 0.7; 0.3 0.4])) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end
ts = 1:the_params.tot_steps÷the_params.save_step

for t in ts
    tadj = t-1
    if tadj % the_params.plot_step != 0
        continue
    end
    local p1 = scatter(Vhistory[t][:, 1], Vhistory[t][:, 2], xlims=[0, 1], ylims=[0, 1], color="#3188CE", legend=:none, marker=:circle, alpha=0.9, markerstrokewidth=0, markersize=2, aspect_ratio=:equal, margin=5Plots.mm, dpi=300)
    scatter!(p1, Phistory[t][:, 1], Phistory[t][:, 2], color=:red, marker=:xcross, alpha=1, markerstrokewidth=7, markersize=10, grid=false)
    savefig(p1, "sims/$(foldername)/t=$(@sprintf("%.3f",tadj/the_params.tot_steps*the_params.tot_time)).png")
end

orders = Dorder.(Vhistory[1:the_params.plot_step:end], Phistory[1:the_params.plot_step:end])
orders = hcat(orders...)

Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(Vhistory[1:the_params.plot_step:end])

p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")


p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders, xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)

if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end
savefig(p2, "sims/$(foldername)/2d_Dhat.png")
savefig(p3, "sims/$(foldername)/2d_Qvv.png")



open("sims/$(foldername)/test.txt", "w") do file
    write(file, "Figure3.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 4 and 6a")
#region Figure 4 and 6a
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=500, # Number of voters
    Np=4, # Number of parties
    Rvv=0.05,# The threshold for a voter consider anothers' opinion
    μvv=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.1, # The threshold for a voter consider anothers' opinion
    μpv=0.8,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0.02, # The standard deviation of the noise added to the voters
    Rvp=0.4, # The threshold for a voter consider anothers' opinion
    μvp=0.02, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.05, # The threshold for a voter consider anothers' opinion
    μpp=0.02, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.002, # The standard deviation of the noise added to the voters
    tot_time=300, # Total time to run the model
    Δt=0.03, # The time step
    plot_step=10, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="periodic"
)


foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.2 0.35 0.6 0.78], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.

(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the data from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:2:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.9, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)

orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)



if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end


savefig(p1, "sims/$(foldername)/convergencenoisy.png")
savefig(p2, "sims/$(foldername)/convergencenoisy_Dhat.png")
savefig(p3, "sims/$(foldername)/convergencenoisy_Qvv.png")


open("sims/$(foldername)/test.txt", "w") do file
    write(file, "Figure5_convergence.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 7")
#region Figure 7
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=500, # Number of voters
    Np=4, # Number of parties
    Rvv=0.05,# The threshold for a voter consider anothers' opinion
    μvv=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.1, # The threshold for a voter consider anothers' opinion
    μpv=0.8,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0.1, # The standard deviation of the noise added to the voters
    Rvp=0.4, # The threshold for a voter consider anothers' opinion
    μvp=0.02, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.05, # The threshold for a voter consider anothers' opinion
    μpp=0.02, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.002, # The standard deviation of the noise added to the voters
    tot_time=300, # Total time to run the model
    Δt=0.03, # The time step
    plot_step=10, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="periodic")


foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.2 0.35 0.6 0.78], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the data from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:2:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.9, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)

orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)



if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end


savefig(p1, "sims/$(foldername)/diffusenoisy.png")
savefig(p2, "sims/$(foldername)/diffusenoisy_Dhat.png")
savefig(p3, "sims/$(foldername)/diffusenoisy_Qvv.png")

open("sims/$(foldername)/test.txt", "w") do file
    write(file, "figure5_divergence.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 8b")
#region Figure 8b
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=200, # Number of voters
    Np=4, # Number of parties
    Rvv=0.6,# The threshold for a voter consider anothers' opinion
    μvv=0.3, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.6, # The threshold for a voter consider anothers' opinion
    μpv=0.3,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0, # The standard deviation of the noise added to the voters
    Rvp=0.6, # The threshold for a voter consider anothers' opinion
    μvp=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.6, # The threshold for a voter consider anothers' opinion
    μpp=0.5, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.0, # The standard deviation of the noise added to the voters
    tot_time=50, # Total time to run the model
    Δt=0.005, # The time step
    plot_step=5, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="reflective")


foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.1, 0.3, 0.5, 0.7], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the data from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:1:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.7, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 250), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)

orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)



if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end


savefig(p1, "sims/$(foldername)/no_convergence_no_noise_evolution.png")
savefig(p2, "sims/$(foldername)/no_convergence_no_noise_evolution_Dhat.png")
savefig(p3, "sims/$(foldername)/no_convergence_no_noise_evolution_Qvv.png")



open("sims/$(foldername)/test.txt", "w") do file
    write(file, "figure6_no_convergence_no_noise.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 8a")
#region Figure 8a
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=200, # Number of voters
    Np=4, # Number of parties
    Rvv=0.6,# The threshold for a voter consider anothers' opinion
    μvv=0.3, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.6, # The threshold for a voter consider anothers' opinion
    μpv=0.3,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0, # The standard deviation of the noise added to the voters
    Rvp=0.6, # The threshold for a voter consider anothers' opinion
    μvp=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.6, # The threshold for a voter consider anothers' opinion
    μpp=0.4, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.0, # The standard deviation of the noise added to the voters
    tot_time=50, # Total time to run the model
    Δt=0.005, # The time step
    plot_step=5, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="reflective"
)


foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.1, 0.3, 0.5, 0.7], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.

(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the dtaa from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:1:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.7, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 250), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)


orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)



if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end


savefig(p1, "sims/$(foldername)/convergence_no_noise_evolution.png")
savefig(p2, "sims/$(foldername)/convergence_no_noise_evolution_Dhat.png")
savefig(p3, "sims/$(foldername)/convergence_no_noise_evolution_Qvv.png")

open("sims/$(foldername)/test.txt", "w") do file
    write(file, "figure6_convergence_no_noise.\n Params:\n")
    write(file, string(the_params))
end
#endregion

println("HANDLING Figure 8c")
#region Figure 8c
the_params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=200, # Number of voters
    Np=4, # Number of parties
    Rvv=0.1,# The threshold for a voter consider anothers' opinion
    μvv=0.3, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.1, # The threshold for a voter consider anothers' opinion
    μpv=0.3,#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=0, # The standard deviation of the noise added to the voters
    Rvp=0.1, # The threshold for a voter consider anothers' opinion
    μvp=0.5, # The rate at which voters change their minds, called the convergence rate
    Rpp=0.1, # The threshold for a voter consider anothers' opinion
    μpp=0.4, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.0, # The standard deviation of the noise added to the voters
    tot_time=50, # Total time to run the model
    Δt=0.005, # The time step
    plot_step=5, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="reflective"
)


foldername = generate_initial_conditions(the_params, (rand(the_params.Nv, 1), reshape([0.1, 0.3, 0.5, 0.7], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
# foldername = generate_initial_conditions(the_params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.

(V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
run_and_save(V0, P0, foldername, the_params) # Run the model with the initial conditions

# Read the data from the foldername
Phistory = load("data/$foldername/results.jld2")["Phistory"]
Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

p1 = scatter(LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Vhistory[1:1:end, 1:the_params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.7, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 250), margin=5Plots.mm, dpi=300)
plot!(p1, LinRange(0.01, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step), Phistory[:, 1:the_params.plot_step:end]', color=:red)


orders = Dorder.(eachcol(Vhistory[:, 1:the_params.plot_step:end]), eachcol(Phistory[:, 1:the_params.plot_step:end]))
orders = hcat(orders...)
Qwanghere(V) = Qwang(V,the_params)
Qwang_orders = Qwanghere.(eachcol(Vhistory[:, 1:the_params.plot_step:end]))
p2 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),orders', label=L"\hat{D}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
plot!(p2, LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),Qwang_orders,label=L"Q_{vv}")

p3 = plot(LinRange(0, the_params.tot_time, the_params.tot_steps ÷ the_params.plot_step),legend=:none,Qwang_orders,label=L"Q_{vv}", xlabel="t", ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)

if !isdir("sims/$foldername")
    mkdir("sims/$foldername")
end


savefig(p1, "sims/$(foldername)/no_convergence_no_noise2_evolution.png")
savefig(p2, "sims/$(foldername)/no_convergence_no_noise2_evolution_Dhat.png")
savefig(p3, "sims/$(foldername)/no_convergence_no_noise2_evolution_Qvv.png")

open("sims/$(foldername)/test.txt", "w") do file
    write(file, "figure6_no_convergence_no_noise2.\n Params:\n")
    write(file, string(the_params))
end
#endregion
