module  plotting_helpers
using LinearAlgebra,Distributions, Distances
export Dorder,Dorder_full, Qwang
    
    function total_order(V,P,params) # Deprecated
        qvv = 0
        for v_i in eachrow(V)
            for v_j in eachrow(V)
                if norm(v_j - v_i) ≤ params.Rvv
                    qvv+=1
                end
            end
        end
        qvv = qvv/(params.Nv^2)
        
        qpp = 0
        for p_i in eachrow(P)
            for p_j in eachrow(P)
                if norm(p_j - p_i) ≤ params.Rpp
                    qpp+=1
                end
            end
        end
        qpp = qpp/(params.Np^2)

        qpv = 0
        for p_α in eachrow(P)
            for v_j in eachrow(V)
                if norm(p_α - v_j) ≤ params.Rpv
                    qpv+=1
                end
            end
        end
        qpv = qpv/(params.Np * params.Nv)

        qvp = 0
        for p_α in eachrow(P)
            for v_j in eachrow(V)
                if norm(p_α - v_j) ≤ params.Rvp
                    qvp+=1
                end
            end
        end
        qvp = qvp/(params.Np * params.Nv)
        
        return qvv,qpp,qpv,qvp
    end
    function C(V)
        # thevar = var(V)
        if length(size(V))==2
            d = size(V)[2]
            thevar = tr(cov(V))
            return min((d/12)/(3*thevar),1)
        else 
            thevar = tr(cov(V))
            return min((1/12)/(3*thevar),1)
        end
    end
    function C(thevar::Number)
        return min((1/12)/(3*thevar),1)
    end

    function C_a(a,V)
        if length(size(V))==2
            d = size(V)[2]
            thevar = tr(cov(V))
            return min((d/12)/(a*thevar),1)
        else 
            thevar = tr(cov(V))
            return min((1/12)/(a*thevar),1)
        end
    end

    function C_a(a,thevar::Number)
        return min((1/12)/(a*thevar),1)
    end

    function Dorder(V,P)
        if ndims(V)==1
            Dvv = pairwise(Euclidean(), V, V) # Voter Voter Distances
            dvv = mean(abs.(Dvv .- round.(Dvv)))

            Dpv = pairwise(Euclidean(), V, P) # Party Voter Distances
            dpv = mean(abs.(Dpv .- round.(Dpv)))

            Dpp = pairwise(Euclidean(), P, P) # Party Party Distances
            dpp = mean(abs.(Dpp .- round.(Dpp)))
            return 1 - 4*(dvv+dpv+dpp)/3
        else
            Dvvs=[]
            Dpvs=[]
            Dpps=[]

            for dim in 1:ndims(V)
                Dvv_dim = pairwise(Euclidean(), V[:,dim], V[:,dim])
                Dpv_dim = pairwise(Euclidean(), V[:,dim], P[:,dim])
                Dpp_dim = pairwise(Euclidean(), P[:,dim], P[:,dim])
                push!(Dvvs,(Dvv_dim .- round.(Dvv_dim)).^2)
                push!(Dpvs,(Dpv_dim .- round.(Dpv_dim)).^2)
                push!(Dpps,(Dpp_dim .- round.(Dpp_dim)).^2)
            end
            dvv = mean(sqrt.(sum(Dvvs)))
            dpv = mean(sqrt.(sum(Dpvs)))
            dpp = mean(sqrt.(sum(Dpps)))
            k = 6 / (log(1+sqrt(2)) + sqrt(2)) # Constant for d2. More complicated for higher d.
            return 1 - k*(dvv+dpv+dpp)/3

        end
    end
    function Dorder_full(V,P)
        if ndims(V)==1
            Dvv = pairwise(Euclidean(), V, V) # Voter Voter Distances
            dvv = mean(abs.(Dvv .- round.(Dvv)))

            Dpv = pairwise(Euclidean(), V, P) # Party Voter Distances
            dpv = mean(abs.(Dpv .- round.(Dpv)))

            Dpp = pairwise(Euclidean(), P, P) # Party Party Distances
            dpp = mean(abs.(Dpp .- round.(Dpp)))
            return @. 1 -  4 * [(dvv+dpv+dpp)/3, dvv, dpv, dpp]
        else
            Dvvs=[]
            Dpvs=[]
            Dpps=[]

            for dim in 1:ndims(V)
                Dvv_dim = pairwise(Euclidean(), V[:,dim], V[:,dim])
                Dpv_dim = pairwise(Euclidean(), V[:,dim], P[:,dim])
                Dpp_dim = pairwise(Euclidean(), P[:,dim], P[:,dim])
                push!(Dvvs,(Dvv_dim .- round.(Dvv_dim)).^2)
                push!(Dpvs,(Dpv_dim .- round.(Dpv_dim)).^2)
                push!(Dpps,(Dpp_dim .- round.(Dpp_dim)).^2)
            end
            dvv = mean(sqrt.(sum(Dvvs)))
            dpv = mean(sqrt.(sum(Dpvs)))
            dpp = mean(sqrt.(sum(Dpps)))
            k = 6 / (log(1+sqrt(2)) + sqrt(2)) # Constant for d2. More complicated for higher d.
            return @. 1 - k* [(dvv+dpv+dpp)/3, dvv, dpv, dpp]

        end
    end
    function Qwang(V,params)
        if ndims(V)==1
            Dvv = pairwise(Euclidean(), V, V) # Voter Voter Distances
            return sum((Dvv .< params.Rvv) .| (Dvv.> (1-params.Rvv)))/params.Nv^2
        else
            Dvv_dims = []
            for dim in 1:ndims(V)
            
                Dvv_dim = pairwise(Euclidean(), V[:,dim], V[:,dim])
                push!(Dvv_dims,(Dvv_dim .- round.(Dvv_dim)).^2) # Applying the periodic boundary condition along each dimension
            end
            return sum(sqrt.(sum(Dvv_dims)) .< params.Rvv) / params.Nv^2 # Apply the Q to each entry.
        end
    end

end