using FreqTables
using StatsBase
using Plots, StatsPlots
using Statistics
using DataFrames
Pkg.build("Conda")
ENV["PYTHON"]=""; Pkg.build("PyCall")
using PyPlot
using RCall
using CSV

include("folder_path/methodsComparison_auxFun.jl")
include("folder_path/method_IFP.jl")
include("folder_path/method_exact.jl")
include("folder_path/method_mcmcsampling.jl")

function compare_approaches(k_range, N_range, percentile_cutoff_range, rho_range, niter)
    results = []
    ii = 1
    iters = length(k_range) * length(N_range) * length(percentile_cutoff_range) * niter * length(rho_range)
    for k in k_range
        for N in N_range
            for percentile_cutoff in percentile_cutoff_range
                for rho in rho_range
                    for iter in 1:niter
                        @info "iteration: $ii of $iters"
                        t = time()
                        ifp = method_IFP(k,k,N,percentile_cutoff, rho, verbose = false)
                        t_ipf = t - time()
                        t = time()
                        exc = method_exact(k,k,N,percentile_cutoff, rho, verbose = false)
                        t_exc = time() - t
                        t = time()
                        try
                            mcm = method_sampling_mcmc(k,k,N,percentile_cutoff, rho, verbose = false)
                            t_mcm = time() - t
                            push!(results,(k, N, percentile_cutoff, rho, ifp, exc, mcm, t_ipf, t_exc, t_mcm))
                        catch
                        end
                        ii += 1
                    end
                end
            end
        end
    end
    results
end

results = compare_approaches(3:6, 1000, 5:5:30, 0.00:0.05:0.50, 100)
results = results[[!isnothing(w) for w in results]]

kk = [w[1] for w in results]
nn = [w[2] for w in results]
pc = [w[3] for w in results]
rr = [w[4] for w in results]
ipf_mae = [w[5].mae for w in results]
ipf_vc = [w[5].vc for w in results]
ipf_bl = [w[5].blanks for w in results]
exa_mae = [w[6].mae for w in results]
exa_vc = [w[6].vc for w in results]
exa_bl = [w[6].blanks for w in results]
mcm_mae = [w[7].mae for w in results]
mcm_vc = [w[7].vc for w in results]
mcm_bl = [w[7].blanks for w in results]
ipf_t = [w[8] for w in results]
exa_t = [w[9] for w in results]
mcm_t = [w[10] for w in results]

results = DataFrame(K = kk, N = nn, PC = pc, RR = rr, IPF_blanks = ipf_bl, IPF_time = ipf_t, IPF_VC = ipf_vc, IPF_MAE = ipf_mae, EXA_blanks = exa_bl, EXA_time = exa_t, EXA_VC = exa_vc, EXA_MAE = exa_mae, MCM_blanks = mcm_bl, MCM_time = mcm_t, MCM_VC = mcm_vc, MCM_MAE = mcm_mae)

CSV.write("folder_path/CellSampling_results.csv", results, delim=";", decimal = ',')
