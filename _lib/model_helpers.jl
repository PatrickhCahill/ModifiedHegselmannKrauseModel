module  model_helpers


################################ Imports ################################
using JLD2
using Dates
using ProgressBars
using Random
using LinearAlgebra

################################ Private Functions ################################
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


################################ Public Functions ################################
export generate_initial_conditions, load_initial_conditions, myrun, run_and_save
function generate_initial_conditions(params) # No intial state so generate a random one
    V0 = rand(params.Nv,params.dimension) # Initial voter beliefs
    P0 = rand(params.Np,params.dimension) # Initial party belcciefs

    filename = "initial_conditions"
    folder = "data/$(string(today()))_run_1"
    if !isdir("$folder")
        mkdir("$folder")
    else
        # Get all the folders in the data directory
        folders = readdir("data")
        # Remove any folders that don't begin with today's date
        folders = [current_folder for current_folder in folders if (length(current_folder) >= length(string(today())))]
    
        folders = [current_folder for current_folder in folders 
                    if current_folder[1:length(string(today()))] == string(today())]
    
        # Only get folders where the beginning is today's date
        folders = [current_folder for current_folder in folders if current_folder[1:length(string(today()))] == string(today())]
        # Get the number after "run_"
        folder_index = [current_folder[length("$(string(today()))_run_1"):end] for current_folder in folders]
        folder_index_int = [0 for _ in folder_index]
        for (index,element) in enumerate(folder_index)
            try
                folder_index_int[index] = parse(Int,element)
            catch e
                println("Warning: Ignoring incorrectly named folder.")
                println(e)
            end
        end
        folder = "data/$(string(today()))_run_$(maximum(folder_index_int)+1)"
        mkdir("$folder")
    end
    # Save the initial conditions
    description = """Initial Conditions generate on $(today()).
        data_file = ":/data/$filename.jld2"
        N = $(params.Nv)
        numP = $(params.Np)
    
        Voters are drawn from a uniform distribution on [0,1]x[0,1].
        Parties are at centres of each quadrant of the unit square. i.e (0.25,0.25), (0.25,0.75), (0.75,0.25), (0.75,0.75).
    
    
    """
    
    
    open("$folder/$filename.txt", "w") do io
        write(io, description)
    end
    
    save("$folder/$filename.jld2", "V0", V0, "P0", P0)

    return  folder[6:end] # Remove the "data/" from the folder name
end

function generate_initial_conditions(params,state::Tuple{Matrix{Float64}, Matrix{Float64}}) # Initial state is given
    (V0,P0) = state
    if size(V0) != (params.Nv,params.dimension)
        throw(ArgumentError("V0 must be $(params.Nv)x$(params.dimension) to match params.Nv"))
    elseif size(P0) != (params.Np,params.dimension)
        throw(ArgumentError("P0 must be $(params.Np)x$(params.dimension) to match params.Np"))
    end

    filename = "initial_conditions"
    folder = "data/$(string(today()))_run_1"
    if !isdir("$folder")
        mkdir("$folder")
    else
        # Get all the folders in the data directory
        folders = readdir("data")
        # Remove any folders that don't begin with today's date
        folders = [current_folder for current_folder in folders if (length(current_folder) >= length(string(today())))]
    
        folders = [current_folder for current_folder in folders 
                    if current_folder[1:length(string(today()))] == string(today())]
    
        # Only get folders where the beginning is today's date
        folders = [current_folder for current_folder in folders if current_folder[1:length(string(today()))] == string(today())]
        # Get the number after "run_"
        folder_index = [current_folder[length("$(string(today()))_run_1"):end] for current_folder in folders]
        folder_index_int = [0 for _ in folder_index]
        for (index,element) in enumerate(folder_index)
            try
                folder_index_int[index] = parse(Int,element)
            catch e
                println("Warning: Ignoring incorrectly named folder.")
                println(e)
            end
        end
        folder = "data/$(string(today()))_run_$(maximum(folder_index_int)+1)"
        mkdir("$folder")
    end
    # Save the initial conditions
    description = """Initial Conditions generate on $(today()).
        data_file = ":/data/$filename.jld2"
        N = $(params.Nv)
        numP = $(params.Np)
    
        Voters are drawn from a uniform distribution on [0,1]x[0,1].
        Parties are at centres of each quadrant of the unit square. i.e (0.25,0.25), (0.25,0.75), (0.75,0.25), (0.75,0.75).
    
    
    """
    
    
    open("$folder/$filename.txt", "w") do io
        write(io, description)
    end
    
    save("$folder/$filename.jld2", "V0", V0, "P0", P0)

    return folder[6:end] # Remove the "data/" from the folder name
end

function generate_initial_conditions(params,state)
    throw(ArgumentError("state must be $(typeof((randn(2,2),randn(2,2))))"))
end

function load_initial_conditions(foldername)
    initialData = load("data/$(foldername)/initial_conditions.jld2")
    V0 = initialData["V0"]
    P0 = initialData["P0"]
    return V0,P0
end

function myrun(V,P,params)
    Vhistory = []
    Phistory = []
    if (lowercase(params.boundary_conditions) ∉ ("reflective","periodic"))
        println("$(params.boundary_conditions) is not reflective or periodic. Assuming periodic")
    end
    # Run the update method until the system converges
    for t = ProgressBar(1:params.tot_steps)
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

function run_and_save(V0,P0,foldername,params)
    V,P, Vhistory,Phistory = myrun(V0,P0,params)
    save("data/$foldername/results.jld2","V",V,"P",P,"Vhistory",Vhistory,"Phistory",Phistory)
end


end # module
