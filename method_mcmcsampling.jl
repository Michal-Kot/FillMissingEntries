using RCall

__jhitandrun__ = RCall.rimport("hitandrun")
using .__jhitandrun__

rgenerate_controls = reval("generate_controls = function(XY, XY_obs, threshold){

  XY_unknown = 1 - sign(XY_obs)

  row_constraints = matrix(0, nrow = nrow(XY_unknown), ncol = prod(dim(XY_obs)))
  for(i in 1:nrow(XY)){
    for(j in 1:ncol(XY)){
      row_constraints[i,((i-1)*ncol(XY)+j)] = XY_unknown[i,j]
    }
  }

  col_constraints = matrix(0, nrow = ncol(XY_unknown), ncol = prod(dim(XY_obs)))
  for(i in 1:nrow(XY)){
    for(j in 1:ncol(XY)){
      col_constraints[j,((i-1)*ncol(XY)+j)] = XY_unknown[i,j]
    }
  }


  A1 = rbind(row_constraints, col_constraints)
   A1 = A1[,apply(A1,2,function(x) sum(x) > 0)]
   A2 = diag(-1, ncol(A1))
   A3 = diag(1, ncol(A1))

   b1 = c(rowSums(XY) - rowSums(XY_obs), colSums(XY) - colSums(XY_obs))
   b2 = rep(-1, nrow(A2))
   b3 = rep(threshold, nrow(A3))

   dirs1 = rep('=', length(b1))
   dirs2 = rep('<=', nrow(A2))
   dirs3 = rep('<=', nrow(A3))

   A = rbind(A1, A2, A3)
   b = c(b1, b2, b3)
   dirs = c(dirs1, dirs2, dirs3)

  to_keep = apply(A,1,sum) != 0

  A = A[to_keep,]
  b = b[to_keep]
  dirs = dirs[to_keep]

  to_keep

  controls = list(constr = A, rhs = as.numeric(b), dir = dirs)

  controls

}")

rmapping = reval("mapping = function(sampling, XY, XY_obs, threshold){
  x = apply(sampling,2,mean)
  XY_unknown = 1 - sign(XY_obs)
  k=1
  for(i in 1:nrow(XY)){
    for(j in 1:ncol(XY)){
      if(XY_unknown[i,j] == 1){
        XY_obs[i,j] = x[k]
        k = k+1
      }
    }
  }

  XY_obs

}")

rminus1 = reval("minus1 = function(ctrls){
  list(constr = ctrls[['constr']][-1,], rhs = ctrls[['rhs']][-1], dir = ctrls[['dir']][-1])
}")

function method_sampling_mcmc(kx, ky, n, percentile_cutoff, rho; verbose=true)
  verbose && @info "generating data"
  data = gendata(kx, ky, n, percentile_cutoff, rho)
  verbose && @info "data has size $(size(data.gen)) with $(count(ismissing, data.obs)) missings"
  A = data.gen
  A_obs = @. ifelse(ismissing(data.obs), 0, data.obs)
  A_obs = complete_trival(A, A_obs, data.threshold)
  if any(A_obs .== 0)
      rXY = robject(A)
      rXY_obs = robject(A_obs)
      rthreshold = robject(data.threshold)
      rctrls = rcall(rgenerate_controls, XY=rXY, XY_obs = rXY_obs, threshold = rthreshold)
      rctrls1 = rcall(rminus1, ctrls = rctrls)
      try
          rsampling = rcall(shakeandbake, constr = rctrls)
          rA_obs  = rcall(rmapping, sampling = rsampling, XY = rXY, XY_obs = rXY_obs, threshold = rthreshold)
          A_obs = rcopy(rA_obs)
          mae = mean_absolute_error(data, A_obs)
          return (data=data, blanks=count(ismissing, data.obs), mae=mae, fit = A_obs, vc = cramer_v(data.gen), full=true)
      catch
          try
              rsampling = rcall(shakeandbake, constr = rctrls, eliminate = false)
              rA_obs  = rcall(rmapping, sampling = rsampling, XY = rXY, XY_obs = rXY_obs, vc = cramer_v(data.gen), threshold = rthreshold)
              A_obs = rcopy(rA_obs)
              mae = mean_absolute_error(data, A_obs)
              return (data=data, blanks=count(ismissing, data.obs), mae=mae, fit = A_obs, vc = cramer_v(data.gen), full=true)
          catch
              rsampling = rcall(shakeandbake, constr = rctrls1, eliminate = false)
              rA_obs  = rcall(rmapping, sampling = rsampling, XY = rXY, XY_obs = rXY_obs, threshold = rthreshold)
              A_obs = rcopy(rA_obs)
              mae = mean_absolute_error(data, A_obs)
              return (data=data, blanks=count(ismissing, data.obs), mae=mae, fit = A_obs, vc = cramer_v(data.gen), full = false)
          end
      end
  else
      mae = mean_absolute_error(data, A_obs)
      return (data=data, blanks=count(ismissing, data.obs), mae=mae, fit = A_obs, vc = cramer_v(data.gen), full=true)
  end
end
