using Distributions

function gendata2(k, n, percentile_cutoff, rho)
    X = rand(1:k, n)
    if rho > 0
        Y = ifelse.(rand(Uniform(0,1),length(X)) .<= rho, X, rand(1:k,length(X)))
    else
        Y = rand(1:k, n)
    end
    A = Matrix(freqtable(X,Y))
    threshold = floor.(Int, percentile(vec(A),percentile_cutoff))
    observed_A = @. ifelse(A > threshold, A, 0)
    independent_A = sum(A,dims = 2) * sum(A,dims = 1) ./ sum(A)
    A = @. ifelse(A==0,1,A)
    (obs=observed_A, gen=A, ind = independent_A, threshold=threshold)
end

function mean_absolute_error(data, prediction)
    sum(abs.(data.gen .- prediction)) / count(ismissing.(data.obs))
end

function complete_trival(a, a_obs, threshold)

    if !all(a .== a_obs)

        a_row_sum = [sum(a[i,:]) for i in 1:size(a)[1]]
        a_col_sum = [sum(a[:,j]) for j in 1:size(a)[2]]

        singleton_row = [size(a)[2] - sum(sign.(a_obs[i,:])) for i in 1:size(a)[1]]
        singleton_col = [size(a)[1] - sum(sign.(a_obs[:,j])) for j in 1:size(a)[2]]

        singleton = any(vcat(singleton_row, singleton_col) .== 1)

        while singleton

        a_row_sum_current = [sum(a_obs[i,:]) for i in 1:size(a)[1]]
        a_col_sum_current = [sum(a_obs[:,j]) for j in 1:size(a)[2]]

            for i in 1:size(a_obs)[1]
                for j in 1:size(a_obs)[2]
                    if singleton_row[i] == 1 && a_obs[i,j] == 0
                        a_obs[i,j] = a_row_sum[i] - a_row_sum_current[i]
                        singleton_row[i] = 0
                        singleton_col[j] -= 1
                        a_row_sum_current = [sum(a_obs[i,:]) for i in 1:size(a)[1]]
                        a_col_sum_current = [sum(a_obs[:,j]) for j in 1:size(a)[2]]
                    end
                    if singleton_col[j] == 1 && a_obs[i,j] == 0
                        a_obs[i,j] = a_col_sum[j] - a_col_sum_current[j]
                        singleton_col[j] = 0
                        singleton_row[i] -= 1
                        a_row_sum_current = [sum(a_obs[i,:]) for i in 1:size(a)[1]]
                        a_col_sum_current = [sum(a_obs[:,j]) for j in 1:size(a)[2]]
                    end
                end
            end

            for i in 1:size(a_obs)[1]
                if (sum(1 .- sign.(a_obs[i,:])) * threshold) == (a_row_sum[i] - a_row_sum_current[i]) && (sum(sign.(a_obs[i,:])) < length(axes(a_obs,2)))
                    for j in 1:size(a_obs)[2]
                        if a_obs[i,j] == 0
                            a_obs[i,j] = threshold
                        end
                    end
                end
            end

            for j in 1:size(a_obs)[2]
                if (sum(1 .- sign.(a_obs[:,j])) * threshold) == (a_col_sum[j] - a_col_sum_current[j]) && (sum(sign.(a_obs[:,j])) < length(axes(a_obs,1)))
                    for i in 1:size(a_obs)[1]
                        if a_obs[i,j] == 0
                            a_obs[i,j] = threshold
                        end
                    end
                end
            end

        singleton_row = [size(a)[2] - sum(sign.(a_obs[i,:])) for i in 1:size(a)[1]]
        singleton_col = [size(a)[1] - sum(sign.(a_obs[:,j])) for j in 1:size(a)[2]]

            singleton = any(vcat(singleton_row, singleton_col) .== 1)

        end
    end

    a_obs

end

function compare_single(N, k, percentile_cutoff, rho = 0, correction = true)
    data_sampled = gendata(rand(2:k),N,percentile_cutoff,rho)
    if !(any(data_sampled.gen .== 0))
        #println((N,k,percentile_cutoff))
        mae_ipf = mean_absolute_error(data_sampled, IFP_with_threshold(data_sampled))
        #println(1)
        mae_mys = mean_absolute_error(data_sampled, mySampling2(data_sampled))
        #println(2)
        mae_knn = mean_absolute_error(data_sampled, knnSampling2(data_sampled))
        #println(3)
        (N, k, percentile_cutoff, mae_ipf, mae_mys, mae_knn, cramer_v(data_sampled.gen), rho)
    end
end

function compare_all(N_range, k_range, percentile_cutoff_range, rho_range, niter, correction)
    results = []
    println(length(N_range) * length(k_range) * length(percentile_cutoff_range) * length(rho_range) * niter)
    for N in N_range
        for k in k_range
            for percentile_cutoff in percentile_cutoff_range
                for rho in rho_range
                    for iter in 1:niter
                        try
                            push!(results, compare_single(N, k, percentile_cutoff, rho, correction))
                        catch e
                        end
                    end
                end
            end
        end
    end
    results
end

function create_correlated(N, k, rho)
    X_k = rand(2:k)
    Y_k = rand(2:k)
    X = rand(1:X_k, N)
    if rho > 0
        Y = generate_correlated(X, rho, Y_k)
    else
        Y = rand(1:Y_k, N)
    end
    A = FreqTables.freqtable(X,Y)
end



function cramer_v(a, corrected = true)
    n = sum(a)
    phi = chi_sq(a) / n
    k = size(a)[2]
    r = size(a)[1]
    if corrected
        phi_hat = max(0, phi - (k-1) * (r-1)/(n-1))
        k_hat = k - (k-1)^2 / (n-1)
        r_hat = r - (r-1)^2 / (n-1)
        sqrt(phi_hat / min(k_hat-1,r_hat-1))
    else
        sqrt(phi / min(k-1,r-1))
    end
end

function chi_sq(a)

    row_sums = [sum(a[i,:]) for i in 1:size(a)[1]]
    col_sums = [sum(a[:,j]) for j in 1:size(a)[2]]

    a_hat = row_sums * transpose(col_sums) ./ sum(a)

    sum((a .- a_hat).^2 ./ a_hat)

end

function tapply(x, index)
    res = Float64[]
    for i in sort(unique(index))
        push!(res, mean(x[(index .== i) .& (.!isnan.(x) )]))
    end
    res
end

function range_tapply(x, index, range)
        res = Float64[]
        ranges = collect(range)
        for r in 1:(length(ranges)-1)
            push!(res,mean(x[(index .>= ranges[r]) .& (index .< ranges[r+1]) .& .!(isnan.(x))]))
        end
        res
end
