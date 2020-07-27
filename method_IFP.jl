function method_IFP(kx, ky, n, percentile_cutoff, rho; verbose=true)

    verbose && @info "generating data"
    data = gendata(kx, ky, n, percentile_cutoff, rho)
    verbose && @info "data has size $(size(data.gen)) with $(count(ismissing, data.obs)) missings"

    A = data.gen
    A_obs = @. ifelse(ismissing(data.obs), 0, data.obs)

    initial_matrix = rand(1:1, size(data.gen)[1], size(data.gen)[2])

    real_row_sums = [sum(data.gen[i,:]) for i in 1:size(data.gen)[1]]
    real_col_sums = [sum(data.gen[:,j]) for j in 1:size(data.gen)[2]]

    fixed_row_sums = [sum((A_obs)[i,:]) for i in 1:size(data.gen)[1]]
    fixed_col_sums = [sum((A_obs)[:,j]) for j in 1:size(data.gen)[2]]

    A_over_threshold = sign.(A_obs)
    A_under_threshold = 1 .- A_over_threshold

    test_max = 100

    while test_max >= 1.0001

        current_row_sums = [sum((initial_matrix .* A_under_threshold)[i,:]) for i in 1:size(data.gen)[1]]
        row_multiplier = (real_row_sums - fixed_row_sums) ./ current_row_sums
        replace!(row_multiplier, NaN => 1)
        initial_matrix = initial_matrix .* row_multiplier

        current_col_sums = [sum((initial_matrix .* A_under_threshold)[:,j]) for j in 1:size(data.gen)[2]]
        col_multiplier = (real_col_sums - fixed_col_sums) ./ current_col_sums
        replace!(col_multiplier, NaN => 1)
        initial_matrix = transpose(transpose(initial_matrix) .* col_multiplier)

        current_row_sums = [sum((initial_matrix .* A_under_threshold .+ data.gen .* A_over_threshold)[i,:]) for i in 1:size(data.gen)[1]]
        row_multiplier = real_row_sums ./ current_row_sums
        current_col_sums = [sum((initial_matrix .* A_under_threshold .+ data.gen .* A_over_threshold)[:,j]) for j in 1:size(data.gen)[2]]
        col_multiplier = real_col_sums ./ current_col_sums

        test_max = maximum(vcat(row_multiplier, col_multiplier))
    end

    A_obs = initial_matrix .* A_under_threshold .+ A_obs

    mae = mean_absolute_error(data, A_obs)

    verbose && @info "mae: $mae"

    return(data=data, blanks=count(ismissing, data.obs), mae=mae, vc = cramer_v(data.gen), fit = A_obs)

end
