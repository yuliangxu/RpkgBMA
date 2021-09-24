library(BayesGPfit)
library(MASS)
library(Matrix)
library(mvtnorm)
library(cluster) 
library(RSpectra)
set.seed(1008)

GP.simulate.curve.fast.new = function(x,poly_degree,a,b,
                                      center=NULL,scale=NULL,max_range=6){
  
  x = cbind(x)
  d = ncol(x)
  
  if(is.null(center)){
    center = apply(x,2,mean)
  }
  c_grids = t(x) - center
  if(is.null(scale)){
    max_grids =pmax(apply(c_grids,1,max),-apply(c_grids,1,min))
    scale=as.numeric(max_grids/max_range)
  }
  
  work_x = GP.std.grids(x,center=center,scale=scale,max_range=max_range)
  Xmat = GP.eigen.funcs.fast(grids=work_x,
                             poly_degree =poly_degree,
                             a =a ,b=b)
  lambda = GP.eigen.value(poly_degree=poly_degree,a=a,b=b,d=d)
  betacoef = rnorm(ncol(Xmat),mean=0,sd=sqrt(lambda))
  f = Xmat%*%betacoef
  return(list(f=f,x=x,work_x=work_x, eigen.func = Xmat, eigen.value = lambda))
}
STGP_generate = function(true.image.a,  true.image.b, 
                         sd.noise.a,   sd.noise.b, grids,
                         n.sample){
  # set true params here:
  n_C = 5
  ca = cb = 1
  zetay = zetam = sample(-2:2,n_C)
  C = matrix(rnorm(n_C * n.sample), ncol = n.sample)
  gamma = 1
  X = rnorm(n.sample)
  
  P = length(true.image.a)
  # eta = rnorm(P)
  # Generate covariance matrix
  a = 0.01;b=1
  norm2 = function(x){sum(abs(x)^1.99)}
  Sigma_k = matrix(NA, nrow=P, ncol=P)
  for(i in 1:P){
    for(j in 1:P){
      Sigma_k[i,j] = exp(-a*(norm2(grids[i,])+norm2(grids[j,])) - b*norm2(grids[i,]-grids[j,]) )
    }
  }
  # eta = cumsum(rnorm(P,0,sqrt(1/P)))
  eta = mvrnorm(n=1, mu = rep(0,P), Sigma = Sigma_k)
  M = ca + true.image.a%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta%*%t(rep(1,n.sample)) + matrix(rnorm(n.sample * P, sd = sd.noise.a), ncol = n.sample)
  Y = cb + t(true.image.b%*%M/sqrt(P)) + gamma*X + t(C)%*%zetay + rnorm(n.sample,mean=0, sd = sd.noise.b)
  if(sum(true.image.a!=0)>n.sample | sum(true.image.b!=0)>n.sample){
    stop("nonzero pixel number exceeds sample size.")
  }
  M.nonzero = t(M)[,true.image.b!=0]/sqrt(P)
  snratio = 1-sum( (Y - cb- gamma*X - t(C)%*%zetay - M.nonzero%*%solve(crossprod(M.nonzero),t(M.nonzero)%*%Y) )^2 )/sum((Y-mean(Y))^2)
  
  
  # M-regression(check alpha)
  M.reduced = M - ca - rep(1,P)%*% t(zetam)%*%C - eta%*%t(rep(1,n.sample))
  ind = which(true.image.a!=0)
  R2.alpha = rep(NA, P)
  for(j in ind){
    M_j = M.reduced[j,]
    R2.alpha[j] = 1-sum( (M_j - X*(sum(X*M_j))/sum(X^2) )^2 )/sum((M_j-mean(M_j))^2)
  }
  
  
  true.ll = NULL
  true.ll$beta = -sum((Y-(cb + t(true.image.b%*%M/sqrt(P)) + gamma*X + t(C)%*%zetay))^2)/2/sd.noise.b^2
  true.ll$alpha = apply(M-(ca + true.image.a%*%t(X) + rep(1,P)%*% t(zetam)%*%C),1,function(x){-sum(x^2)/2/sd.noise.a^2})
  return(list(Y = Y, M = M, ca = ca, cb=cb, beta=true.image.b, alpha = true.image.a, X=X,C=C,
              zetay = zetay, zetam = zetam,gamma = gamma,eta=eta,
              sigma_y = sd.noise.b, sigma_m = sd.noise.a,
              true.ll = true.ll, snratio = snratio, R2.alpha = na.omit(R2.alpha)))
}

STGP_generate_theta = function(true.image.a,  true.image.b, 
                         sd.noise.a,   sd.noise.b, grids,Q,lambda,
                         n.sample){
  # set true params here:
  n_C = 5
  ca = cb = 1
  zetay = zetam = sample(-2:2,n_C)
  C = matrix(rnorm(n_C * n.sample), ncol = n.sample)
  gamma = 1
  X = rnorm(n.sample)
  
  P = length(true.image.a)
  L = dim(Q)[2]
  # eta = rnorm(P)
  # Generate covariance matrix
  a = 0.01;b=1
  # norm2 = function(x){sum(abs(x)^1.99)}
  # Sigma_k = matrix(NA, nrow=P, ncol=P)
  # for(i in 1:P){
  #   for(j in 1:P){
  #     Sigma_k[i,j] = exp(-a*(norm2(grids[i,])+norm2(grids[j,])) - b*norm2(grids[i,]-grids[j,]) )
  #   }
  # }
  # eta = cumsum(rnorm(P,0,sqrt(1/P)))
  theta_eta = matrix(rnorm(n*L),nrow=L) #L by n
  # Q_t = big_transpose(Q)
  eta = Q%*%theta_eta
  true_theta_alpha = t(true.image.a%*%Q)
  true_theta_beta = t(true.image.b%*%Q)
  # L = dim(Q)[2]
  # true_theta_alpha = rep(0,L)
  # true_theta_beta = rep(0,L)
  M = ca + Q%*%true_theta_alpha%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta + matrix(rnorm(n.sample * P, sd = sd.noise.a), 
                                                                                                      ncol = n.sample)
  Y = cb + t(M/sqrt(P))%*%ST_fun(Q%*%true_theta_beta,lambda,50,F) + gamma*X + t(C)%*%zetay + rnorm(n.sample,mean=0, sd = sd.noise.b)
  # if(sum(true.image.a!=0)>n.sample | sum(true.image.b!=0)>n.sample){
  #   stop("nonzero pixel number exceeds sample size.")
  # }
  # M.nonzero = t(M)[,true.image.b!=0]/sqrt(P)
  # snratio = 1-sum( (Y - cb- gamma*X - t(C)%*%zetay - M.nonzero%*%solve(crossprod(M.nonzero),t(M.nonzero)%*%Y) )^2 )/sum((Y-mean(Y))^2)
  # 
  
  # # M-regression(check alpha)
  # M.reduced = M - ca - rep(1,P)%*% t(zetam)%*%C - eta%*%t(rep(1,n.sample))
  # ind = which(true.image.a!=0)
  # R2.alpha = rep(NA, P)
  # for(j in ind){
  #   M_j = M.reduced[j,]
  #   R2.alpha[j] = 1-sum( (M_j - X*(sum(X*M_j))/sum(X^2) )^2 )/sum((M_j-mean(M_j))^2)
  # }
  
  
  true.ll = NULL
  true.ll$beta = -sum((Y-(cb + t(M/sqrt(P))%*%ST_fun(Q%*%true_theta_beta,lambda,50,F) + gamma*X + t(C)%*%zetay))^2)/2/sd.noise.b^2
  true.ll$alpha = apply(M-(ca + true.image.a%*%t(X) + rep(1,P)%*% t(zetam)%*%C),1,function(x){-sum(x^2)/2/sd.noise.a^2})
  return(list(Y = Y, M = M, ca = ca, cb=cb, beta=true.image.b, alpha = true.image.a, X=X,C=C,
              zetay = zetay, zetam = zetam,gamma = gamma,theta_eta=theta_eta,
              theta_alpha = true_theta_alpha, theta_beta = true_theta_beta,
              sigma_y = sd.noise.b, sigma_m = sd.noise.a,
              true.ll = true.ll ))
              # snratio = snratio, R2.alpha = na.omit(R2.alpha)))
}

STGP_generate_theta_block = function(true.image.a,  true.image.b, 
                               sd.noise.a,   sd.noise.b, grids,Q,lambda,region_idx,L_all,
                               n.sample){
  # set true params here:
  n_C = 2
  ca = 0; cb = 0
  zetay = runif(n_C,-2,2); zetam = runif(n_C,-2,2)/1e2
  C = matrix(rnorm(n_C * n.sample), ncol = n.sample)*5
  gamma = 0.2
  X = rnorm(n.sample)*5
  true.image.b = true.image.b*2
  
  P = length(true.image.a)
  L = dim(Q[[1]])[2]
  num_block = length(Q)
  
  theta_eta = vector(mode = "list", length = num_block)
  eta = matrix(NA, nrow=p,ncol=n.sample)
  beta_test = rep(0,p)
  alpha_test = rep(0,p)
  true_theta_alpha = rep(NA, sum(L_all))
  true_theta_beta = rep(NA, sum(L_all))
  sigma_eta=0.1
  theta_eta_final = NULL
  for(m in 1:num_block){
    print(paste("m=",m))
    L = L_all[m]
    L_end = cumsum(L_all)[m]
    theta_eta[[m]] = matrix(rnorm(n.sample*L)*sigma_eta,nrow=L, ncol = n.sample) #L by n
    idx = region_idx[[m]]
    print(paste("dim(Q[[m]])=",dim(Q[[m]]),";dim(theta_eta[[m]])=",dim(theta_eta[[m]])))
    eta[idx,] = Q[[m]]%*%theta_eta[[m]] # p_block by n
    theta_eta_final = rbind(theta_eta_final,theta_eta[[m]])
    true_theta_alpha[(L_end-L+1):L_end] = t(true.image.a[idx]%*%Q[[m]])
    true_theta_beta[(L_end-L+1):L_end] = t(true.image.b[idx]%*%Q[[m]])
    beta_test[idx] = Q[[m]]%*%true_theta_beta[(L_end-L+1):L_end]
    alpha_test[idx] = Q[[m]]%*%true_theta_alpha[(L_end-L+1):L_end]
  }
  beta_test_ST = (beta_test - sign(beta_test)*lambda)*(abs(beta_test)>lambda)
  alpha_test_ST = (alpha_test - sign(alpha_test)*lambda)*(abs(alpha_test)>lambda)
  M = ca + alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta + matrix(rnorm(n.sample * P, sd = sd.noise.a), 
                                                                                 ncol = n.sample)
  Y = cb + t(M/sqrt(P))%*% beta_test_ST + gamma*X + t(C)%*%(zetay) + rnorm(n.sample,mean=0, sd = sd.noise.b)
  
  # R^2
  # R^2 for M at each location: higher R^2 better accuracy
  snratio = NULL
  M.nonzero = t(M)[,beta_test_ST!=0]/sqrt(P)
  snratio$R2.beta = 1-sum( (Y - cb- gamma*X - t(C)%*%zetay - M.nonzero%*%solve(crossprod(M.nonzero),t(M.nonzero)%*%Y) )^2 )/sum((Y-mean(Y))^2)
  # resY = Y-(cb + t(M/sqrt(P))%*% beta_test_ST + gamma*X + t(C)%*%(zetay))
  snratio$sgn.beta = sd( t(M/sqrt(P))%*% beta_test_ST )/sd( Y-(cb  + gamma*X + t(C)%*%(zetay)) )
  snratio$gamma = sd(gamma*X)/sd( Y-(cb + t(M/sqrt(P))%*% beta_test_ST + t(C)%*%(zetay)) )
  snratio$zetay = sd(t(C)%*%zetay)/sd( Y-(cb + t(M/sqrt(P))%*% beta_test_ST + gamma*X )  )
  snratio$Ymodel = sd( t(M/sqrt(P))%*% beta_test_ST + gamma*X )/sd(Y)
   

  # M-regression(check alpha)
  M.reduced = M - (ca + alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta)
  signal = alpha_test_ST%*%t(X)
  snratio$sgn.alpha = rep(NA, P)
  for(j in 1:P){
    M_j = M.reduced[j,]
    f_j = signal[j,]
    snratio$sgn.alpha[j] = sd(f_j)/sd(M_j)
  }
  
  true.ll = NULL
  true.ll$beta = -sum((Y-(cb + t(M/sqrt(P))%*% beta_test_ST + gamma*X + t(C)%*%zetay))^2)/2/sd.noise.b^2
  logll_M = 0
  for(m in 1:num_block){
    Q_t = t(Q[[m]])
    idx = region_idx[[m]]
    p_i = length(idx)
    L = L_all[m]
    L_end = cumsum(L_all)[m]
    L_idx = (L_end-L+1):L_end
    M_star_m = Q_t%*%(M[idx,] - as.matrix(rep(1,p_i))%*%t(t(C)%*%zetam)-alpha_test_ST[idx]%*%t(X)) - theta_eta_final[L_idx,];
    logll_M = logll_M -0.5*(norm(M_star_m,"f"))^2/sd.noise.a^2
  }
  
  norm_M = norm( M -( alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C)-eta, type="f")
  true.ll$alpha = logll_M
  return(list(Y = Y, M = M, ca = ca, cb=cb, beta=true.image.b, alpha = true.image.a, X=X,C=C,
              zetay = zetay, zetam = zetam,gamma = gamma,theta_eta=theta_eta_final,
              theta_alpha = true_theta_alpha, theta_beta = true_theta_beta,
              sigma_y = sd.noise.b, sigma_m = sd.noise.a,snratio=snratio,
              true.ll = true.ll, beta_test_ST=beta_test_ST,alpha_test_ST=alpha_test_ST,sigma_eta=sigma_eta ))
}

FDR = function(active_region, true_region){
  sum(active_region!=0 & true_region==0)/sum(active_region!=0)
}
Precision = function(active_region, true_region){
  mean(I(active_region!=0) == I(true_region!=0))
}
Power = function(active_region, true_region){
  sum(active_region !=0 & true_region!=0)/sum(true_region!=0)
}

matern_kernel = function(x,y,nu,l=1){
  d = sqrt(sum((x-y)^2))/l
  y = 2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d)^nu*besselK(sqrt(2*nu)*d,nu)
  return(y)
}
generate_approx_matern_basis = function(grids, region_idx_list,scale = 2,nu = 2){
  # compute centers
  num_block = length(region_idx_list)
  clusters = vector("list",num_block)
  all_centers = vector("list",num_block)
  pt = sapply(region_idx_list, length)
  Lt = sapply(pt,function(x){ceiling(x*0.028)})
  grids_std = NULL
  for(i in 1:num_block){
    print(paste("region/block=",i))
    # standardize grids
    # equal spaced grid centers
    grids_c = sweep(grids[region_idx_list[[i]],],2,apply(grids[region_idx_list[[i]],],2,mean),"-")
    grids_scale = pmax(apply(grids_c,1,max),-apply(grids_c,1,min))
    max_range = sd(grids_scale)
    grids_scale = grids_scale/max_range
    grids_std_i = sweep(grids_c,1,grids_scale ,"/")
    # grids_std = rbind(grids_std,grids_std_i)
    grids_std = rbind(grids_std,grids_c)
    # get centers, equally spaced centers
    clusters[[i]] = kmeans(grids_std_i,centers = Lt[i],iter.max = 20)
    kmfit.centers = fitted(clusters[[i]],"centers")
    kmfit.labels = fitted(clusters[[i]],"classes")
    centers = kmfit.centers[!duplicated(kmfit.labels),]
    all_centers[[i]] = centers
  }
  
  kernel_centers = vector("list",num_block)
  kernel_eg = vector("list",num_block)
  kernel_grids = vector("list",num_block)
  temp = vector("list",num_block)
  Phi = vector("list",num_block)
  Phi_D = vector("list",num_block)
  Phi_Q = vector("list",num_block)
  for(i in 1:num_block){
    print(paste("i=",i))
    kernel = matrix(NA,nrow = Lt[i], ncol=Lt[i])
    k_grid = matrix(NA,nrow = pt[i], ncol=Lt[i])
    for(l in 1:Lt[i]){
      kernel[l,] = apply(all_centers[[i]],1,matern_kernel,y=all_centers[[i]][l,],nu = nu,l=scale)
      k_grid[,l] = apply(grids_std[region_idx_list[[i]],],1,matern_kernel,y=all_centers[[i]][l,],nu = nu,l=scale)
    }
    diag(kernel) = 1
    if(any(c(anyNA(kernel),anyNA(k_grid)))){
      print("any(c(anyNA(kernel),anyNA(k_grid))) = TRUE")
      break
    }
    kernel_centers[[i]] = kernel
    kernel_grids[[i]] = k_grid
    eg = eigen(kernel)
    kernel_eg[[i]] = eg
    temp[[i]] = k_grid%*%eg$vectors
    Phi[[i]] = sqrt(Lt[i])*sweep(temp[[i]],2,eg$values,"/")
    if(anyNA(Phi[[i]])){
      print("anyNA(Phi[[i]]) = True")
      break
    }
    Phi_D[[i]] = eg$values
    # Phi_Q[[i]] = Phi[[i]]
    Phi_Q[[i]] = qr.Q(qr(Phi[[i]]))
  }
  return(list(kernel_centers =kernel_centers,
              all_centers=all_centers,
                   kernel_grids = kernel_grids,
                   kernel_eg = kernel_eg,
                   Phi_D = Phi_D,
                   region_idx_block = region_idx_list,
                   Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
}

generate_matern_basis2 = function(grids, region_idx_list,scale = 2,nu = 1/5,L = 231){
  num_block = length(region_idx_list)
  Phi_D = vector("list",num_block)
  Phi_Q = vector("list",num_block)
  Lt = NULL; pt = NULL
  for(i in 1:num_block){
    p_i = length(region_idx_list[[i]])
    kernel_mat = matrix(NA,nrow = p_i, ncol=p_i)
    for(l in 1:p_i){
      kernel_mat[l,] = apply(grids[region_idx_list[[i]],],1,matern_kernel,y=grids[region_idx_list[[i]],][l,],nu = nu,l=scale)
    }
    diag(kernel_mat) = 1
    K = eigs_sym(kernel_mat,L)
    K_QR = qr(K$vectors)
    Phi_Q[[i]] = qr.Q(K_QR )
    Phi_D[[i]] = K$values
    Lt = c(Lt, length(Phi_D[[i]]))
    pt = c(pt, dim(Phi_Q[[i]])[1])
  }
  return(list(Phi_D = Phi_D,
              region_idx_block = region_idx_list,
              Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
}

generate_sq_basis = function(grids, region_idx_list,a = 0.01,b = 10){
  num_block = length(region_idx_list)
  Phi_D = vector("list",num_block)
  Phi_Q = vector("list",num_block)
  Lt = NULL; pt = NULL
  for(i in 1:num_block){
    GP = GP.simulate.curve.fast.new(x=grids[region_idx_list[[i]],], a=0.01 ,b=10,poly_degree=20) # try to tune b, increase for better FDR
    K_esq = GP$eigen.func
    K_QR = qr(K_esq)
    Phi_Q[[i]] = qr.Q(K_QR)
    Phi_D[[i]] = GP$eigen.value
    Lt = c(Lt, length(Phi_D[[i]]))
    pt = c(pt, dim(Phi_Q[[i]])[1])
  }
  return(list(Phi_D = Phi_D,
              region_idx_block = region_idx_list,
              Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
}

# data2 = STGP_generate(1,image2,sd.noise = 1,n.sample=100)
# data1$snratio
# range(data1$R2.alpha)
# GP.plot.curve(list(f=data1$beta,x=grids), main = "")
# GP.plot.curve(list(f=data1$alpha,x=grids), main = "")



# =======test========== #
if(0){
  data1 = STGP_generate(true.image.a = image1.a,
                        true.image.b = image1.b,
                        sd.noise.a = 0.1, sd.noise.b = 0.1,
                        grids = grids,
                        n.sample=1000)
  data1$snratio
  
  data2 = STGP_generate_theta(true.image.a = image1.a,
                              true.image.b = image1.b,
                              sd.noise.a = 0.1, sd.noise.b = 0.1,
                              grids = grids,
                              n.sample=1000)
  test = matrix(rep(0,J),nrow = p_n)
  j = i = 1:p_n
  for(i in 1:p_n){
    for(j in 1:p_n){
      if( (i-p_n/4-2)^2/9+(j-p_n/2)^2/4 <=2 | (i-3*p_n/4)^2/4+(j-p_n/2)^2/9 <=2 ){
        test[i,j]=1
      }
    }
  }
  test.image = as.vector(test)
  GP.plot.curve(list(f=test.image, x=grids), main = "test image")
  data1 = STGP_generate(true.image.a = test.image,
                        true.image.b = image1.b,
                        sd.noise.a = 0.1, sd.noise.b = 1,
                        n.sample=100)
}
