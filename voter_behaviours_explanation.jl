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

    boundary_conditions::String # Defines which boundary conditions are to be used
end

#region Figure2_Voter Behaviours
runset = [
    [0.1,0.1,1,0.03,0.2,  [0.2, 0.53, 0.86]], # Political Base
    [0.1,0.35,1,0.03,0.2, [0.2, 0.53, 0.86]], # Swing voters only
    [0.1,0.1,1,0.01,0.2,  [0.2, 0.53, 0.86]], # Competition by noise reduction for swing voters
    [0.1,0.1,1,0.03,0.2,  [0.4, 0.53, 0.76]], # Disaffected voters
    [0.3,0.1,1,0.03,0.2,  [0.2, 0.53, 0.86]], # Revolution
    [0.1,0.23,1,0.03,0.2,  [0.2, 0.53, 0.86]], # Competition by R_pv 


]
for set in runset
    params = Params(
        dimension=1, # The dimension of the model | Should be 1 or 2
        Nv=1000, # Number of voters
        Np=3, # Number of parties
        Rvv=set[1],# The threshold for a voter consider anothers' opinion
        μvv=1, # The rate at which voters change their minds, called the convergence rate
        Rpv=set[2], # The threshold for a voter consider anothers' opinion
        μpv=set[3],#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
        σv=set[4], # The standard deviation of the noise added to the voters
        Rvp=set[5], # The threshold for a voter consider anothers' opinion
        μvp=0.03, # The rate at which voters change their minds, called the convergence rate
        Rpp=0.1, # The threshold for a voter consider anothers' opinion
        μpp=0.01, # The rate at which voters change their minds, calledx the convergence rate
        σp=0.005, # The standard deviation of the noise added to the voters
        tot_time=100, # Total time to run the model
        Δt=0.02, # The time step
        plot_step=20, # The period between which the plots are produced
        save_step=1, # How often to save the model
        ensembleInitialNoise=0.02, # The initial noise of the ensemble
        poll_noise_std=0.04, # The variance of the noise in the poll
        n_members=100,
        inflation_factor=1,
        poll_freq=10,
        boundary_conditions="periodic"
    )
    println(params)

    foldername = generate_initial_conditions(params, (rand(params.Nv, 1), reshape(set[6], (:, 1)))) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.
    # foldername = generate_initial_conditions(params) # Generate the initial conditions with random voters and parties and save it to a specific location. We may parse through a specific initial state if we want.

    (V0, P0) = load_initial_conditions(foldername) # Load the initial conditions from the foldername
    run_and_save(V0, P0, foldername, params) # Run the model with the initial conditions

    # Read the data from the foldername
    Phistory = load("data/$foldername/results.jld2")["Phistory"]
    Vhistory = load("data/$foldername/results.jld2")["Vhistory"]
    Phistory = hcat(Phistory...)
    Vhistory = hcat(Vhistory...)

    p1 = scatter(LinRange(0, params.tot_time, params.tot_steps ÷ params.plot_step), Vhistory[1:2:end, 1:params.plot_step:end]', color="#3188CE", legend=:none, xlabel="t", marker=:circle, alpha=0.9, markerstrokewidth=0, markersize=0.5, ylim=[0, 1], size=(600, 300), margin=5Plots.mm, dpi=300)
    plot!(p1, LinRange(0.01, params.tot_time, params.tot_steps ÷ params.plot_step), Phistory[:, 1:params.plot_step:end]', color=:red)







    if !isdir("sims/$foldername")
        mkdir("sims/$foldername")
    end

    savefig(p1, "sims/$(foldername)/noisyconvergence.png")




    open("sims/$(foldername)/test.txt", "w") do file
        write(file, "Figure2.\n Params:\n")
        write(file, string(params))
    end
end
#endregion


