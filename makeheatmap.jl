using Parameters
using Statistics
using LinearAlgebra
using HypothesisTests
using Random
using Distances
Random.seed!(3)


function hegselmann_krause_update(V,P,params)
    dV = zeros(size(V))
    dP = zeros(size(P))

    sampleN_v = params.Nv ÷ 10
    samplejs = Vector(1:params.Nv)
    V_i = zeros(2)
    P_α = zeros(2)
    δV = zeros(2)
    δP = zeros(2)
    
    @views for voter_index in 1:params.Nv
        shuffle!(samplejs)

        @. V_i = V[voter_index,:]
        dV[voter_index,1] += randn()*params.σv*sqrt(params.Δt) # Add noise to the voter
        dV[voter_index,2] += randn()*params.σv*sqrt(params.Δt) # Add noise to the voter

        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - V_i
            @. δV = δV - round(δV)
            if 0 < norm(δV) < params.Rvv
                @. dV[voter_index, :] += δV * params.μvv/sampleN_v*params.Δt
            end
        end
        
        # Attempt to avoid eachrow loop
        # diffs_i = (V[samplejs,:]' .- V_i)'
        # dV[voter_index,:] += params.μvv/sampleN_v * sum(mapslices(x->φ(x/params.Rvv) * x, diffs_i, dims=1)) # Update the voter
        
        for i in 1:params.Np
            @. δP = P[i, :] - V_i
            @. δP = δP - round(δP)

            if 0 < norm(δP) < params.Rpv
                @. dV[voter_index, :] += δP * params.μpv/params.Np*params.Δt
            end
        end
    end

    @views for party_index in 1:params.Np
        shuffle!(samplejs)

        @. P_α = P[party_index,:]
        dP[party_index,1] = randn()*params.σp*sqrt(params.Δt) # Add noise to the party
        dP[party_index,2] = randn()*params.σp*sqrt(params.Δt) # Add noise to the party
        
        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - P_α
            @. δV = δV - round(δV)

            if 0 < norm(δV) < params.Rvp
                @. dP[party_index,:] += δV * params.μvp/sampleN_v*params.Δt
            end
        end
        
        for i in 1:params.Np
            @. δP = P[i, :] - P_α
            @. δP = δP - round(δP)

            if 0 < norm(δP) < params.Rpp
                @. dP[party_index,:] -= params.μpp/params.Np * δP*params.Δt
            end
        end
    end


    V = mod.(V .+ dV, 1)
    P = mod.(P .+ dP, 1)

    return V,P
end

function hk1d(V,P,params)
    dV = zeros(size(V))
    dP = zeros(size(P))

    sampleN_v = params.Nv ÷ 10
    samplejs = Vector(1:params.Nv)
    V_i = zeros(1)
    P_α = zeros(1)
    δV = zeros(1)
    δP = zeros(1)
    
    @views for voter_index in 1:params.Nv
        shuffle!(samplejs)

        @. V_i = V[voter_index,:]
        dV[voter_index,1] += randn()*params.σv*sqrt(params.Δt) # Add noise to the voter

        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - V_i
            @. δV = δV - round(δV)
            if 0 < norm(δV) < params.Rvv
                @. dV[voter_index, :] += δV * params.μvv/sampleN_v*params.Δt
            end
        end
        
        # Attempt to avoid eachrow loop
        # diffs_i = (V[samplejs,:]' .- V_i)'
        # dV[voter_index,:] += params.μvv/sampleN_v * sum(mapslices(x->φ(x/params.Rvv) * x, diffs_i, dims=1)) # Update the voter
        
        for i in 1:params.Np
            @. δP = P[i, :] - V_i
            @. δP = δP - round(δP)

            if 0 < norm(δP) < params.Rpv
                @. dV[voter_index, :] += δP * params.μpv/params.Np*params.Δt
            end
        end
    end

    @views for party_index in 1:params.Np
        shuffle!(samplejs)

        @. P_α = P[party_index,:]
        dP[party_index,1] = randn()*params.σp*sqrt(params.Δt) # Add noise to the party
        
        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - P_α
            @. δV = δV - round(δV)

            if 0 < norm(δV) < params.Rvp
                @. dP[party_index,:] += δV * params.μvp/sampleN_v*params.Δt
            end
        end
        
        for i in 1:params.Np
            @. δP = P[i, :] - P_α
            @. δP = δP - round(δP)

            if 0 < norm(δP) < params.Rpp
                @. dP[party_index,:] -= params.μpp/params.Np * δP*params.Δt
            end
        end
    end

    V = mod.(V .+ dV, 1)
    P = mod.(P .+ dP, 1)
    return V,P
end

function hk1dreflective(V,P,params)
    dV = zeros(size(V))
    dP = zeros(size(P))

    sampleN_v = params.Nv ÷ 10
    samplejs = Vector(1:params.Nv)
    V_i = zeros(1)
    P_α = zeros(1)
    δV = zeros(1)
    δP = zeros(1)
    
    @views for voter_index in 1:params.Nv
        shuffle!(samplejs)

        @. V_i = V[voter_index,:]
        dV[voter_index,1] += randn()*params.σv*sqrt(params.Δt) # Add noise to the voter

        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - V_i
            if 0 < norm(δV) < params.Rvv
                @. dV[voter_index, :] += δV * params.μvv/sampleN_v*params.Δt
            end
        end
        
        # Attempt to avoid eachrow loop
        # diffs_i = (V[samplejs,:]' .- V_i)'
        # dV[voter_index,:] += params.μvv/sampleN_v * sum(mapslices(x->φ(x/params.Rvv) * x, diffs_i, dims=1)) # Update the voter
        
        for i in 1:params.Np
            @. δP = P[i, :] - V_i

            if 0 < norm(δP) < params.Rpv
                @. dV[voter_index, :] += δP * params.μpv/params.Np*params.Δt
            end
        end
    end

    @views for party_index in 1:params.Np
        shuffle!(samplejs)

        @. P_α = P[party_index,:]
        dP[party_index,1] = randn()*params.σp*sqrt(params.Δt) # Add noise to the party
        
        for s in samplejs[1:sampleN_v] 
            @. δV = V[s, :] - P_α

            if 0 < norm(δV) < params.Rvp
                @. dP[party_index,:] += δV * params.μvp/sampleN_v*params.Δt
            end
        end
        
        for i in 1:params.Np
            @. δP = P[i, :] - P_α

            if 0 < norm(δP) < params.Rpp
                @. dP[party_index,:] -= params.μpp/params.Np * δP*params.Δt
            end
        end
    end
    function reflective(v, dv)
        if v+dv <= 1 && v+dv >= 0
            return v+dv
        elseif v+dv > 1
            2- (v+dv)
        else
            return -(v+dv)
        end
    end


    V = reflective.(V,dV)
    P = reflective.(P,dP)
    return V,P
end

function Dorder(V,P)
    Dvv = pairwise(Euclidean(), V, V) # Voter Voter Distances
    dvv = mean(abs.(Dvv .- round.(Dvv)))

    Dpv = pairwise(Euclidean(), V, P) # Party Voter Distances
    dpv = mean(abs.(Dpv .- round.(Dpv)))

    Dpp = pairwise(Euclidean(), P, P) # Party Party Distances
    dpp = mean(abs.(Dpp .- round.(Dpp)))

    return [(dvv+dpv+dpp)/3,dvv,dpv,dpp]
end

function check_convergence(vector)
    n = length(vector)
    slice_size = div(n, 5)
    adf_result = nothing
    # Step through 10% slices
    vals = []
    for i in 1:5
        current_slice = vector[(i-1)*slice_size + 1:(i * slice_size)]

        if length(unique(current_slice))==1
            push!(vals,current_slice[1])
        else
            try
                adf_result = ADFTest(current_slice, :constant,1)
                if pvalue(adf_result) < 0.05 # Typically, a p-value < 0.05 indicates rejection of null hypothesis (converged)
                    push!(vals,mean(current_slice))
                end
            catch e
                # Open "nohup.out" in append mode and write the error message
                open("nohup.out", "a") do file
                    println(file, "A convergence error was caught: ", e)
                    println(file, params)
                end
            end
        end
    end
    if length(vals)>=1
        return (vals[end],true)
    end
    # If none of the slices have converged, return the last 10% slice
    last_slice = vector[(end - slice_size + 1):end]
    return (mean(last_slice), false, pvalue(adf_result))
end

function myrun(V,P,params)
    Vhistory = []
    Phistory = []
    if (lowercase(params.boundary_conditions) ∉ ("reflective","periodic"))
        println("$(params.boundary_conditions) is not reflective or periodic. Assuming periodic")
    end
    # Run the update method until the system converges
    for t = (1:params.tot_steps)
        # Add the current state to the history
        if (t % params.save_step == 0)
            push!(Vhistory, copy(V))
            push!(Phistory, copy(P))
        end


        if lowercase(params.boundary_conditions)=="reflective"
            V,P = hk1dreflective(V,P,params)
        elseif params.dimension == 1
            V,P = hk1d(V,P,params)
        elseif params.dimension == 2
            V,P = hegselmann_krause_update(V,P,params)
        end
    end
    return V,P, Vhistory, Phistory
end 


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
params = Params(
    dimension=1, # The dimension of the model | Should be 1 or 2
    Nv=2000, # Number of voters
    Np=3, # Number of parties
    Rvv=parse(Float64, ARGS[1]),# The threshold for a voter consider anothers' opinion
    μvv=1, # The rate at which voters change their minds, called the convergence rate
    Rpv=0.5, # The threshold for a voter consider anothers' opinion
    μpv=parse(Float64, ARGS[2]),#parse(Float64, ARGS[1]), # The rate at which voters change their minds, called the convergence rate
    σv=parse(Float64, ARGS[3]), # The standard deviation of the noise added to the voters
    Rvp=0.3, # The threshold for a voter consider anothers' opinion
    μvp=0., # The rate at which voters change their minds, called the convergence rate
    Rpp=0.05, # The threshold for a voter consider anothers' opinion
    μpp=0.0, # The rate at which voters change their minds, calledx the convergence rate
    σp=0.0, # The standard deviation of the noise added to the voters
    tot_time=500, # Total time to run the model
    Δt=0.1, # The time step
    plot_step=10, # The period between which the plots are produced
    save_step=1, # How often to save the model
    ensembleInitialNoise=0.02, # The initial noise of the ensemble
    poll_noise_std=0.04, # The variance of the noise in the poll
    n_members=100,
    inflation_factor=1,
    poll_freq=10,
    boundary_conditions="periodic"
)

V = rand(params.Nv, 1)
P = reshape([0.2, 0.53, 0.86], (:, 1))
V,P, Vhistory, Phistory = myrun(V,P,params)
Phistory = hcat(Phistory...)
Vhistory = hcat(Vhistory...)

new_orders = hcat(Dorder.(eachcol(Vhistory)[1:5:end],eachcol(Phistory)[1:5:end])...)
conv_result = check_convergence(new_orders[1,:])

println("$(params.Rvv), $(params.μvv), $(params.μpv), $(params.σv), $(conv_result[1]), $(conv_result[2])")
