const PATH_TO_MINIZINC = "C:/Program Files/MiniZinc/minizinc.exe"

function gendata(kx, ky, n, percentile_cutoff, rho=0)
    X = rand(1:kx, n)
    Y = rand(1:kx, n)
    Y = ifelse.(rand(Uniform(0,1), length(X)) .< rho, X, Y)
    A = Matrix(freqtable(X,Y))
    threshold = floor.(Int, percentile(vec(A),percentile_cutoff))
    independent_A = sum(A, dims=2) * sum(A, dims=1) / sum(A)
    observed_A = @. ifelse(A > threshold, A, missing)
    (obs=observed_A, gen=A, threshold=threshold, ind = independent_A)
end

function listoptions(data)
    locs = findall(ismissing, data.obs)
    open("listallsolutions.mzn", "w") do io
        for loc in locs
            println(io, "var 0..$(data.threshold): a$(loc[1])_$(loc[2]);")
        end
        for i in axes(data.obs, 1)
            print(io, "constraint ")
            l = []
            for j in axes(data.obs, 2)
                if ismissing(data.obs[i,j])
                    push!(l, "a$(i)_$j")
                else
                    push!(l, data.obs[i, j])
                end
            end
            println(io, join(l, "+"), "=", sum(data.gen[i, :]),";")
        end
        for j in axes(data.obs, 2)
            print(io, "constraint ")
            l = []
            for i in axes(data.obs, 1)
                if ismissing(data.obs[i,j])
                    push!(l, "a$(i)_$j")
                else
                    push!(l, data.obs[i, j])
                end
            end
            println(io, join(l, "+"), "=", sum(data.gen[:, j]),";")
        end
        println(io, "solve satisfy;")
        print(io, "output [")
        print(io, join(["show(a$(loc[1])_$(loc[2]))" for loc in locs], ",\" \","))
        println(io, "];")
    end
    df = DataFrame([Int[] for _ in 1:length(locs)],
         [Symbol("a", loc[1], "_", loc[2]) for loc in locs])
    res = read(`$PATH_TO_MINIZINC -a listallsolutions.mzn`, String)
    res_lines = split(res, "\n")
    for i in 1:2:length(res_lines)-2
        push!(df, parse.(Int, split(res_lines[i])))
    end
    return df
end



function compare(data, df)
    agg_df = by(stack(df, :), :variable, :value => mean)
    #closest_idx = argmin(sum(abs.(transpose(convert(Array, df)) .- agg_df.value_mean), dims = 2))[1]
    #agg_df.exact_case = vec(convert(Array,df[closest_idx,:]))
    #println("Closest to the mean is $closest_idx")
    agg_df.actual = [data.gen[parse.(Int, split(String(idx)[2:end], "_"))...] for idx in agg_df.variable]
    return agg_df
end

mse_compare(agg_df) = mean(@. (agg_df.value_mean - agg_df.actual) ^ 2)
mae_compare(agg_df) = mean(@. abs(agg_df.value_mean - agg_df.actual))


function method_exact(k, n, percentile_cutoff,rho; verbose=true)
    verbose && @info "generating data"
    data = gendata(k, k, n, percentile_cutoff,rho)
    verbose && @info "data has size $(size(data.gen)) with $(count(ismissing, data.obs)) missings"
    df = listoptions(data)
    verbose && @info "possible assignments: $(nrow(df))"
    agg_df = compare(data, df)
    agg_df.difference = agg_df.actual - agg_df.value_mean
    mae = mae_compare(agg_df)
    verbose && @info "mae: $mae"
    return (data=data, blanks=count(ismissing, data.obs), options=df, mae=mae, flatten_df = agg_df, vc = cramer_v(data.gen), fit = agg_df.value_mean)
end

using Distributions
using FreqTables
using StatsBase
using DataFrames

exact_variance = [[method_exact(rand(5:6),1000,m,rand(Uniform(0,0.5)), verbose=false).flatten_df.difference for x in 1:100] for m in 10:10:50]


#function run_experiments(kx, ky, n, percentile_cutoff, iters)
#    df = DataFrame(blanks=Int[], options=Int[], mse=Float64[])
#    for i in 1:iters
#        @info "iteration: $i of $iters"
#        _, blanks, options, mse = one_experiment(kx, ky, n, percentile_cutoff, verbose=false)
#        push!(df, (blanks, nrow(options), mse))
#    end
#    df
#end

#@time res = run_experiments(7, 7, 1000, 30, 500);
#res_agg1 = by(res, :blanks, :mse => length, :mse => mean, sort=true)
#plot(res_agg1.blanks, res_agg1.mse_mean)
#res.options_g = cut(res.options, 10)
#res_agg2 = by(res, :options_g, :mse => length, :mse => mean, sort=true)
#plot(res_agg2.options_g, res_agg2.mse_mean)
