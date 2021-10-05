
## Data generation ---------------------------------------------------------


simulate_data<-function(case,n_obs=3000) {

  length_tot<-500
  grid_s<-grid_t<-seq(0,1,length.out = length_tot)



  # generate X --------------------------------------------------------------

  domain<-c(0,1)
  n_basis_x<-32     #random chosen between 10 and 50
  X_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_x)
  X_coef<-matrix(rnorm(n_obs*n_basis_x),nrow = n_basis_x,ncol = n_obs )
  X_fd<-fd(X_coef,X_basis)
  X<-eval.fd(grid_s,X_fd)
  X<- as.matrix(t(as.matrix(rproc2fdata(n_obs,t=grid_s+1,sigma="brownian",par.list=list("scale"=1))[[1]])))
  # Generate ERROR ----------------------------------------------------------


  n_basis_eps<-20 #random chosen between 10 and 50
  eps_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_eps)
  eps_coef<-matrix(rnorm(n_obs*n_basis_eps),nrow = n_basis_eps,ncol = n_obs )
  eps_fd<-fd(eps_coef,eps_basis)
  Eps<-t(eval.fd(grid_t,eps_fd))


  # Define beta -----------------------------------------------------------

  if(case=="Scenario I"){
    cat("Scenario I")
    beta<-function(s,t){

      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')

      z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
      z_matrix}
      else{
        z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
        z_matrix}
    }



  }
  if(case=="Scenario II"){
    cat("Scenario II")
    beta<-function(s,t){
      a=0.25
      b=0.25
      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')
      z<- -(((d$s-0.5)/a)^2 + ((d$t-0.5)/b)^2) +1
      z[z<0]<-0
      z_matrix<-matrix(z,nrow=length(s))
      z_matrix}
      else{
        z<- -(((s)/a)^2 + ((d)/b)^2) + 1
        z[z<0]<-0
        z}
    }}
  if(case=="Scenario III"){
    cat("Scenario III")

    beta<-function(s,t){
      a<-0.05
      b<-0.5
      f_1<-function(s,t){b*(1-s)*sin(10*pi*(s-a-1+sqrt(1-(t-0.5)^2)))}
      f_2<-function(s,t){b*sin(10*pi*(s+a+1-sqrt(1-(t-0.5)^2)))}

      z<-matrix(0,length_tot,length_tot)
      for (ii in 1:length(grid_t)) {
        t<-grid_t[ii]
        s_0_1<-grid_s[grid_s>( a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_0_2<-grid_s[grid_s<( -a+sqrt(1-(t-0.5)^2))&grid_s>=0.5]
        s_n0_1<-grid_s[grid_s<=(a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_n0_2<-grid_s[grid_s>=(-a+sqrt(1-(t-0.5)^2))&grid_s>0.5]
        z_i<-c(f_1(s_n0_1,t),rep(0,length(s_0_1)),rep(0,length(s_0_2)),f_2(s_n0_2,t))
        z[ii,]=z_i
      }
      return(t(z))
    }
  }

  if(case=="Scenario IV"){
    cat("Scenario IV")
    beta<-function(s,t){
      a=0.5
      b=0.5
      c=0.5
      d=0.5
      f_1<-function(s,t){((t-0.5)/c)^3+((s-0.5)/d)^3+((t-0.5)/b)^2 + ((s-0.5)/a)^2+5}
      z<- outer(s,t,f_1)
      z
    }
  }
  if(case=="Historical"){
    cat("Historical")
    beta<-function(grid_s,grid_t){

      f_1<-function(s,t){0.5*sin(30*s)+0.5*cos(30*t)}
      grid<-expand.grid(s=grid_s,t=grid_t)
      z<-f_1(grid$s,grid$t)
      ind_0<-which(grid$s>grid$t)
      z[ind_0]<-0
      matrix(z,length(grid_s),length(grid_t))

    }
  }
  if(case=="Concurrent"){
    cat("Concurrent")
    beta<-function(grid_s,grid_t){

      f_1<-function(s,t){0.5*sin(30*s)+0.5*cos(30*t)}
      grid<-expand.grid(s=grid_s,t=grid_t)
      z<-f_1(grid$s,grid$t)
      ind_0<-which(grid$s!=grid$t)
      z[ind_0]<-0
      matrix(z,length(grid_s),length(grid_t))

    }
  }

  if(case=="Concurrent"){G<-t(eval.basis(grid_s,X_basis))%*%beta(grid_s,grid_t)}
  else{
    G<-(1/length(grid_s))*t(eval.basis(grid_s,X_basis))%*%beta(grid_s,grid_t)
  }
  Y_parz<-t(X_coef)%*%G
   Y_parz<-(1/length(grid_s))*t(X)%*%beta(grid_s,grid_t)
  signal_to_noise_ratio<-40
  ##for each observation SN=(sigma^2_signal/sigma^2_error)
  if(case=="Scenario I"){Y = Y_parz + Eps%*%diag(colVars(Eps)^(-1/2))}
  else{
    k <- sqrt((colVars(Y_parz)+max(colVars(Y_parz)))/(signal_to_noise_ratio*colVars(Eps)))

    Y = Y_parz + Eps%*%diag(k)
  }

  out<-list(X=X,
            Y=t(Y),
            X_fd=X_fd,
            Eps=Eps,
            beta_matrix=beta(grid_s,grid_t))

  return(out)
}




## B SLASSO----------------------------------
slasso.fr.cv<-function(Y_fd,X_fd,basis_s,basis_t,K=10,ss_rule_par=0,lam_opt_method="min",
                       lambdas_L=seq(-2,1,by=0.25),lambdas_s=seq(-4,-1),lambdas_t=seq(-4,1),B0=NA,max_iterations=2000){
  
  n_obs<-dim(X_fd$coefs)[2]
  W_X<-eval.penalty(basis_s)
  W_Y<-eval.penalty(basis_t)
  R_X<-eval.penalty(basis_s,2)
  R_Y<-eval.penalty(basis_t,2)
  domain_s<-basis_s$rangeval
  domain_t<-basis_t$rangeval
  n_basis_s<-basis_s$nbasis
  n_basis_t<-basis_t$nbasis
  orders<-n_basis_s-length(basis_s$params)
  ordert<-n_basis_t-length(basis_t$params)
  breaks_s<-basis_s$params
  breaks_t<-basis_t$params
  ext_break_s<-c(rep(domain_s[1],orders),breaks_s,rep(domain_s[2],orders))
  ext_break_t<-c(rep(domain_s[1],ordert),breaks_t,rep(domain_s[2],ordert))
  weights_s<-diff(ext_break_s,lag=orders)/orders
  weights_t<-diff(ext_break_t,lag=ordert)/ordert
  weights_mat<-weights_s%o%weights_t
  weights_vec<-vec(weights_mat)
  ##Centering data
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  X_fd_cen<-center.fd(X_fd)
  Y_fd_cen<-center.fd(Y_fd)
  
  X_new<-inprod(X_fd_cen,basis_s)
  Y_new<-inprod(Y_fd_cen,basis_t)
  
  env <- new.env()
  #Matrix initialization
  env[["n_basis_x"]] <- n_basis_s
  env[["n_basis_y"]] <- n_basis_t
  env[["Y_newc"]] <- Y_new
  env[["X_newc"]] <- X_new
  env[["W_X"]] <-W_X
  env[["W_Y"]] <-W_Y
  env[["R_X"]] <-R_X
  env[["R_Y"]] <-R_Y
  
  if(is.na(B0)[1]){
    B_basis<-ginv(t(X_new)%*%X_new)%*%t(X_new)%*%Y_new%*%solve(W_Y)
  }
  else{
    B_basis<-B0
  }
  
  
  
  
  
  log_lambda_L_vec<-lambdas_L
  log_lambda_s_vec<-lambdas_s
  log_lambda_t_vec<-lambdas_t
  
  lambda_L_vec<-10^log_lambda_L_vec
  lambda_s_vec<-10^log_lambda_s_vec
  lambda_t_vec<-10^log_lambda_t_vec
  
  comb_list<-expand.grid(lambda_L_vec,lambda_s_vec,lambda_t_vec)
  
  parr_func<-function(ii){
    
    lambdas<-comb_list[ii,]
    lambda_L<-as.numeric(lambdas[1])
    lambda_s<-as.numeric(lambdas[2])
    lambda_t<-as.numeric(lambdas[3])
    
    env[["lambda_x_opt"]] <-lambda_s
    env[["lambda_y_opt"]] <-lambda_t
    ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
    split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
    inpr_vec<-list()
    
    
    for (ll in 1:K) {
      #cat(ll)
      Y_i<-Y_fd_cen[split_vec[[ll]]]
      X_i<-X_new[split_vec[[ll]],]
      Y_minus<-Y_new[-split_vec[[ll]],]
      X_minus<-X_new[-split_vec[[ll]],]
      env[["Y_newc"]] <- Y_minus
      env[["X_newc"]] <- X_minus
      
      output <- lbfgs_2(objective(), gradient(), B_basis, lambda = lambda_L,weights = weights_vec,environment=env,max_iterations = max_iterations, invisible = 1)
      B_par<-matrix(output$par,n_basis_s,n_basis_t)
      Y_hat<-fd(t(X_i%*%B_par),basis_t)
      inpr_vec[[ll]]<-mean(diag(inprod(Y_i-Y_hat,Y_i-Y_hat)))
      
    }
    env[["Y_newc"]] <- Y_new
    env[["X_newc"]] <- X_new
    fit_full <- lbfgs_2(objective(), gradient(), B_basis, lambda = lambda_L,weights = weights_vec,environment=env,max_iterations =100, invisible = 1)
    B_full<-matrix(fit_full$par,nrow=n_basis_s,ncol = n_basis_t)
    Beta_full<-bifd(B_full,basis_s,basis_t)
    per_0<-get_per_0(Beta_full)
    mean<-mean(unlist(inpr_vec))
    sd<-sd(unlist(inpr_vec))/sqrt(K)
    out<-as.numeric(list(mean=mean,
                         sd=sd,
                         per_0=per_0))
    return(out)
    
    
  }
  
  
  cores <- detectCores()
  vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_func,mc.cores = cores)
  par<-sapply(vec_par,"[[",1)
  sd<-sapply(vec_par,"[[",2)
  per_0<-sapply(vec_par,"[[",3)
  
  ##Lambdas optimum
  
  lambdas_opt<-as.numeric(comb_list[max(which(par<=min(par))),])
  
  if(lam_opt_method!="min"){
    
    
    len_s<-length(lambdas_s)
    len_t<-length(lambdas_t)
    len_L<-length(lambdas_L)
    
    lambda_opt_i<-matrix(0,length(ss_rule_par),3)
    for(kkk in 1:length(ss_rule_par)){
      ###Lambda L
      mat_CV<-matrix(0,len_s*len_t,len_L)
      mat_CV_sd<-matrix(0,len_s*len_t,len_L)
      for(ii in 1:len_L){
        mat_CV[,ii]<-par[which(comb_list[,1]==lambda_L_vec[ii])]
        mat_CV_sd[,ii]<-sd[which(comb_list[,1]==lambda_L_vec[ii])]
      }
      if(lam_opt_method=="sd_mean"){
        means_CV_L<-colMeans(mat_CV)
        sd_CV_L<-colMeans(mat_CV_sd)
        ind_opt_L<-max(which(abs(means_CV_L-min(means_CV_L))<=ss_rule_par[kkk]*sd[max(which(par==min(means_CV_L)))]))
        lambda_opt_L[kkk]<-lambda_L_vec[ind_opt_L] 
      }
      else if(lam_opt_method=="sd_min"){
        min_CV_L<-colMins(mat_CV)
        mat_st<-list()
        for(ii in 1:len_L){
          mat_st[[ii]]<-comb_list[which(par==min_CV_L[ii]),]
        }
        sd_CV_L<-colMins(mat_CV_sd)
        ind_opt_L<-max(which(abs(min_CV_L-min(min_CV_L))<=ss_rule_par[kkk]*sd[max(which(par==min(min_CV_L)))]))
        lambda_opt_i[kkk,]<-as.numeric(mat_st[[ind_opt_L]])
      }
      
      
    }
    lambdas_opt<-rbind(lambdas_opt,lambda_opt_i)
  }
  
  out<-list(vec_par=vec_par,
            lambdas_opt=lambdas_opt,
            CV=par,
            CV_sd=sd,
            per_0=per_0,
            comb_list=comb_list,
            X=X_fd,
            Y=Y_fd,
            type="SLASSO")
  return(out)
  
}

slasso.fr<-function(Y_fd,X_fd,basis_s,basis_t,
                    lambdas_L,lambdas_s,lambdas_t,B0=NA,...){
  
  n_obs<-dim(X_fd$coefs)[2]
  W_X<-eval.penalty(basis_s)
  W_Y<-eval.penalty(basis_t)
  R_X<-eval.penalty(basis_s,2)
  R_Y<-eval.penalty(basis_t,2)
  domain_s<-basis_s$rangeval
  domain_t<-basis_t$rangeval
  n_basis_s<-basis_s$nbasis
  n_basis_t<-basis_t$nbasis
  orders<-n_basis_s-length(basis_s$params)
  ordert<-n_basis_t-length(basis_t$params)
  breaks_s<-basis_s$params
  breaks_t<-basis_t$params
  ext_break_s<-c(rep(domain_s[1],orders),breaks_s,rep(domain_s[2],orders))
  ext_break_t<-c(rep(domain_s[1],ordert),breaks_t,rep(domain_s[2],ordert))
  weights_s<-diff(ext_break_s,lag=orders)/orders
  weights_t<-diff(ext_break_t,lag=ordert)/ordert
  weights_mat<-weights_s%o%weights_t
  weights_vec<-vec(weights_mat)
  ##Centering data
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  X_fd_cen<-center.fd(X_fd)
  Y_fd_cen<-center.fd(Y_fd)
  
  X_new<-inprod(X_fd_cen,basis_s)
  Y_new<-inprod(Y_fd_cen,basis_t)
  
  env <- new.env()
  #Matrix initialization
  env[["n_basis_x"]] <- n_basis_s
  env[["n_basis_y"]] <- n_basis_t
  env[["Y_newc"]] <- Y_new
  env[["X_newc"]] <- X_new
  env[["W_X"]] <-W_X
  env[["W_Y"]] <-W_Y
  env[["R_X"]] <-R_X
  env[["R_Y"]] <-R_Y
  
  if(is.na(B0)){
    B_basis<-ginv(t(X_new)%*%X_new)%*%t(X_new)%*%Y_new%*%solve(W_Y)
  }
  else{
    B_basis<-B0
  }
  
  
  
  
  lambda_L<-10^lambdas_L
  env[["lambda_x_opt"]] <- 10^lambdas_s
  env[["lambda_y_opt"]] <-10^lambdas_t
  cat("SLASSO:",c(lambda_L, 10^lambdas_s, 10^lambdas_t),"     ")
  output <- lbfgs_2(objective(), gradient(), B_basis, environment=env,lambda = lambda_L,weights = weights_vec,...)
  B_par<-matrix(output$par,nrow=n_basis_s,ncol = n_basis_t)
  Beta_hat_fd<-bifd(B_par,basis_s,basis_t)
  
  #Intercept
  X_mean_new<-inprod(X_mean,basis_s)
  alpha<-Y_mean-fd(t(X_mean_new%*%B_par),basis_t)
  
  out<-list(B=B_par,
            Beta_hat_fd=Beta_hat_fd,
            alpha=alpha,
            lambdas_L=lambdas_L,
            lambdas_s=lambdas_s,
            lambdas_t=lambdas_t,
            X=X_fd,
            Y=Y_fd,
            type="SLASSO")
  return(out)
  
}

# 
# 

# 
# 
# ## B spline n basis regularization ramesy ----------------------------------
# 
# 
# fr.tru.cv<-function(Y_fd,X_fd,K=10,nbasiss=seq(4,10,by=1),nbasist=seq(4,10,by=1),basis_type="bspline"){
#   
#   X_fd_cen=center.fd(X_fd)
#   Y_fd_cen=center.fd(Y_fd)
#   
#   domainx<-X_fd_cen$basis$rangeval
#   domainy<-Y_fd_cen$basis$rangeval
#   if(basis_type=="bspline"){
#     create_fun=create.bspline.basis
#   }
#   else if(basis_type=="fourier"){
#     create_fun=create.fourier.basis
#   }
#   
#   n_basis_x_seq<-nbasiss
#   n_basis_y_seq<-nbasist
#   n_obs<-length(X_fd_cen$coefs[1,])
#   
#   comb_list<-expand.grid(n_basis_x_seq,n_basis_y_seq)
#   inpr_vec<-list()
#   parr_tru<-function(ii){
#     
#     basis<-as.numeric(comb_list[ii,])
#     n_basis_x<-basis[1]
#     n_basis_y<-basis[2]
#     basis_x <- create_fun(domainx,nbasis = n_basis_x)
#     basis_y <- create_fun(domainy,nbasis=n_basis_y)
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     
#     for(ll in 1:K){
#       
#       Y_i<-Y_fd_cen[split_vec[[ll]]]
#       X_i<-X_fd_cen[split_vec[[ll]],]
#       X_minus<-X_fd_cen[-split_vec[[ll]],]
#       Y_minus<-Y_fd_cen[-split_vec[[ll]],]
#       mod<-fregre.basis.fr(X_minus,Y_minus,basis.s =basis_x,basis.t = basis_y, lambda.s = 0,lambda.t = 0)
#       Y_hat<-predict(mod,X_i)
#       inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
#     }
#     mean<-mean(unlist(inpr_vec))
#     sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
#     out<- list(mean=mean,
#                sd=sd)
#     return(out)
#     
#   }
#   
#   cores <- detectCores()
#   
#   vec_par<-lapply(seq(1,length(comb_list[,1])),parr_tru)#,mc.cores = cores)
#   par<-sapply(vec_par,"[[",1)
#   sd<-sapply(vec_par,"[[",2)
#   
#   
#   
#   n_basis_opt<-as.numeric(comb_list[max(which(par==min(par,na.rm = T))),])
#   
#   
#   
#   out<-list(par_opt=n_basis_opt,
#             CV=par,
#             CV_sd=sd,
#             comb_list=comb_list,
#             X=X_fd,
#             Y=Y_fd,
#             type="TRU")
#   
#   return(out)
#   
# }
# 
# fr.tru<-function(Y_fd,X_fd,nbasiss=5,nbasist=5,basis_type="bspline",...){
#   
#   n_obs<-length(X_fd$coefs[1,])
#   
#   if(basis_type=="bspline"){
#     create_fun=create.bspline.basis
#   }
#   else if(basis_type=="fourier"){
#     create_fun=create.fourier.basis
#   }
#   
#   domainx<-X_fd$basis$rangeval
#   domainy<-Y_fd$basis$rangeval
#   
#   basis_x <- create_fun(domainx,nbasis = nbasiss)
#   basis_y <- create_fun(domainy,nbasis=nbasist)
#   
#   mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, lambda.s = 0,lambda.t = 0)
#   
#   B<-mod$beta.estbifd$coefs
#   Beta_hat_fd<-mod$beta.estbifd
#   cat("TRU:",nbasiss,nbasist,"     ")
#   
#   out<-list(B=B,
#             Beta_hat_fd=Beta_hat_fd,
#             nbasiss=nbasiss,
#             nbasist=nbasist,
#             mod=mod,
#             type="TRU")
#   
#   return(out)
#   
# }
# 
# 
# ## B splines penalized regression Ramsey 2005 -------------------------------
# 
# 
# fr.usc.cv<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=10^seq(5,15,by=1),lambdas_t=10^seq(5,15,by=1)){
#   
#   n_obs<-length(X_fd$coefs[1,])
#   
#   
#   inpr_vec<-list()
#   comb_list<-expand.grid(lambdas_s,lambdas_t)
#   
#   
#   
#   parr_smooth<-function(ii){
#     
#     lambdas<-as.numeric(comb_list[ii,])
#     lambda_s<-lambdas[1]
#     lambda_t<-lambdas[2]
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     
#     for(ll in 1:K){
#       
#       Y_i<-Y_fd[split_vec[[ll]]]
#       X_i<-X_fd[split_vec[[ll]],]
#       X_minus<-X_fd[-split_vec[[ll]],]
#       Y_minus<-Y_fd[-split_vec[[ll]],]
#       mod<-fregre.basis.fr(X_minus,Y_minus,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s,lambda.t = lambda_t)
#       Y_hat<-predict(mod,X_i)
#       inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
#     }
#     mean<-mean(unlist(inpr_vec))
#     sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
#     out<- list(mena=mean,
#                sd=sd)
#     
#     return(out)
#     
#   }
#   cores <- detectCores()
#   vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_smooth,mc.cores = cores)
#   par<-sapply(vec_par,"[[",1)
#   sd<-sapply(vec_par,"[[",2)
#   
#   l_opt<-as.numeric(comb_list[max(which(par==min(par))),])
#   
#   
#   lambda_s_opt<-l_opt[1]
#   lambda_t_opt<-l_opt[2]
#   mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s_opt,lambda.t = lambda_t_opt)
#   
#   B<-mod$beta.estbifd$coefs
#   Beta_hat_fd<-mod$beta.estbifd
#   
#   out<-list(B=B,
#             Beta_hat_fd=Beta_hat_fd,
#             lambda_s_opt=lambda_s_opt,
#             lambda_t_opt=lambda_t_opt,
#             CV=par,
#             CV_sd=sd,
#             comb_list=comb_list,
#             type="SMOOTH")
#   return(out)
# }
# fr.usc<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=0,lambdas_t=0){
#   
#   n_obs<-length(X_fd$coefs[1,])
#   mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambdas_s,lambda.t = lambdas_t)
#   
#   B<-mod$beta.estbifd$coefs
#   Beta_hat_fd<-mod$beta.estbifd
#   l_opt<-c(lambdas_s,lambdas_t)
#   cat("SMOOTH:",l_opt,"     ")
#   out<-list(B=B,
#             Beta_hat_fd=Beta_hat_fd,
#             lambda_s=lambdas_s,
#             lambda_t=lambdas_t,
#             mod=mod,
#             type="SMOOTH")
#   return(out)
# }
# 
# 
# 
# ## PCA regression YAO ----------------------------------
# fr.PCA.cv<-function(Y_fd,X_fd,K=10){
#   
#   X_fd_cen=center.fd(X_fd)
#   Y_fd_cen=center.fd(Y_fd)
#   X_mean<-mean.fd(X_fd)
#   Y_mean<-mean.fd(Y_fd)
#   n_obs<-length(X_fd$coefs[1,])
#   domain_s=X_fd$basis$rangeval
#   domain_t=Y_fd$basis$rangeval
#   length_grid=500
#   grid_s<-seq(domain_s[1],domain_s[2],length.out = length_grid)
#   grid_t<-seq(domain_t[1],domain_t[2],length.out = length_grid)
#   delta_t<-((domain_t[2]-domain_t[1])/length_grid)
#   delta_s<-((domain_s[2]-domain_s[1])/length_grid)
#   X_fd_eval<-t(eval.fd(grid_s,X_fd_cen))
#   Y_fd_eval<-t(eval.fd(grid_t,Y_fd_cen))
#   basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
#   basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
#   grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
#   grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
#   X_list<-split(t(X_fd_eval),rep(1:ncol(t(X_fd_eval)), each = length_grid))
#   Y_list<-split(t(Y_fd_eval),rep(1:ncol(t(Y_fd_eval)), each = length_grid))
#   
#   
#   PCA_x<-FPCA(X_list,grid_list_s,optns = list(useBinnedData="OFF"))
#   PCA_y<-FPCA(Y_list,grid_list_t,optns = list(useBinnedData="OFF"))
#   n_comp_x<-seq(3,PCA_x$selectK-2,by=2)
#   n_comp_y<-seq(3,PCA_y$selectK-2,by=2)
#   
#   comb_list<-expand.grid(n_comp_x,n_comp_y)
#   inpr_vec<-numeric()
#   
#   parr_PCA<-function(ii){
#     
#     n_components<-as.numeric(comb_list[ii,])
#     n_comp_x<-n_components[1]
#     n_comp_y<-n_components[2]
#     
#     
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     CV_ll=list()
#     
#     for(ll in 1:K){
#       
#       Y_i<-Y_fd_cen[split_vec[[ll]]]
#       X_i<-X_fd_eval[split_vec[[ll]],]
#       X_minus<-X_fd_eval[-split_vec[[ll]],]
#       Y_minus<-Y_fd_eval[-split_vec[[ll]],]
#       X_fd_minus<-X_fd_cen[-split_vec[[ll]]]
#       Y_fd_minus<-Y_fd_cen[-split_vec[[ll]]]
#       
#       grid_list_s<-lapply(1:ncol(t(X_minus)), function(ii)grid_s)
#       grid_list_t<-lapply(1:ncol(t(Y_minus)), function(ii)grid_t)
#       
#       X_list<-split(t(X_minus),rep(1:ncol(t(X_minus)), each = nrow(t(X_minus))))
#       Y_list<-split(t(Y_minus),rep(1:ncol(t(Y_minus)), each = nrow(t(Y_minus))))
#       
#       PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
#       PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
#       
#       var_xy<-var.fd(X_fd_minus,Y_fd_minus)
#       var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
#       sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
#       rho<-PCA_x$lambda
#       beta_mat<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
#       Y_hat_mat<-delta_s*X_i%*%beta_mat
#       Y_hat<-fd(t(Y_hat_mat),basis_int_t)
#       CV_ll[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
#     }
#     mean<-mean(unlist(CV_ll))
#     sd<-sd(unlist(CV_ll))/sqrt(n_obs)
#     out<- list(mean=mean,
#                sd=sd)
#     return(out)
#     
#     
#   }
#   
#   cores <- detectCores()
#   vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_PCA,mc.cores = cores)
#   par<-sapply(vec_par,"[[",1)
#   sd<-sapply(vec_par,"[[",2)
#   
#   n_comps_opt<-as.numeric(comb_list[min(which(par==min(par))),])
#   n_comp_x<-n_comps_opt[1]
#   n_comp_y<-n_comps_opt[2]
#   
#   cat("PCA:",n_comps_opt,"     ")
#   grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
#   grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
#   X_list<-split(t(X_fd_eval),rep(1:ncol(t(X_fd_eval)), each = length_grid))
#   Y_list<-split(t(Y_fd_eval),rep(1:ncol(t(Y_fd_eval)), each = length_grid))
#   PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
#   PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
#   
#   var_xy<-var.fd(X_fd_cen,Y_fd_cen)
#   var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
#   sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
#   
#   rho<-PCA_x$lambda
#   
#   beta_mat<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
#   basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
#   basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
#   Beta_fd<-bifd(as.matrix(beta_mat),basis_int_s,basis_int_t)
#   
#   out<-list(Beta_hat_fd=Beta_fd,
#             n_comp_opt_x=n_comp_x,
#             n_comp_opt_y=n_comp_y,
#             CV=par,
#             CV_sd=sd,
#             comb_list=comb_list,
#             X=X_fd,
#             Y=Y_fd,
#             type="PCA")
#   return(out)
#   
# }
# 
# fr.PCA<-function(Y_fd,X_fd,ncomps,ncompt){
#   
#   X_fd_cen=center.fd(X_fd)
#   Y_fd_cen=center.fd(Y_fd)
#   X_mean<-mean.fd(X_fd)
#   Y_mean<-mean.fd(Y_fd)
#   n_obs<-length(X_fd$coefs[1,])
#   domain_s=X_fd$basis$rangeval
#   domain_t=Y_fd$basis$rangeval
#   length_grid=500
#   grid_s<-seq(domain_s[1],domain_s[2],length.out = length_grid)
#   grid_t<-seq(domain_t[1],domain_t[2],length.out = length_grid)
#   delta_t<-((domain_t[2]-domain_t[1])/length_grid)
#   delta_s<-((domain_s[2]-domain_s[1])/length_grid)
#   X_fd_eval<-t(eval.fd(grid_s,X_fd_cen))
#   Y_fd_eval<-t(eval.fd(grid_t,Y_fd_cen))
#   basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
#   basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
#   grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
#   grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
#   
#   n_comps<-c(ncomps,ncompt)
#   n_comp_x<-n_comps[1]
#   n_comp_y<-n_comps[2]
#   
#   cat("PCA:",n_comps,"     ")
#   X_list<-split(t(X_fd_eval),rep(1:n_obs, each = length_grid))
#   Y_list<-split(t(Y_fd_eval),rep(1:n_obs, each = length_grid))
#   PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
#   PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
#   
#   var_xy<-var.fd(X_fd_cen,Y_fd_cen)
#   var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
#   sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
#   
#   rho<-PCA_x$lambda
#   
#   B<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
#   basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
#   basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
#   Beta<-bifd(as.matrix(B),basis_int_s,basis_int_t)
#   #Intercept
#   X_mean_new<-inprod(X_mean,Beta$sbasis)
#   alpha<-Y_mean-fd(t(X_mean_new%*%B),Beta$tbasis)
#   
#   
#   out<-list(Beta_hat_fd=Beta,
#             B=B,
#             alpha=alpha,
#             ncompx=n_comp_x,
#             ncompy=n_comp_y,
#             X=X_fd,
#             Y=Y_fd,
#             type="PCA")
#   return(out)
#   
# }
# 
# 
# ## Ridge vantini canale ----------------------------------------------------
# fr.ridge.cv<-function(Y_fd,X_fd,basiss,basist,K=10,alpha_seq=NA){
#   
#   X_fd_cen=center.fd(X_fd)
#   Y_fd_cen=center.fd(Y_fd)
#   X_mean<-mean.fd(X_fd)
#   Y_mean<-mean.fd(Y_fd)
#   n_obs<-length(X_fd$coefs[1,])
#   domain_s=X_fd$basis$rangeval
#   domain_t=Y_fd$basis$rangeval
#   
#   X_new<-inprod(X_fd_cen,basiss)
#   Y_new<-inprod(Y_fd_cen,basist)
#   
#   W_s<-eval.penalty(basiss)
#   W_t<-eval.penalty(basist)
#   
#   R_s<-eval.penalty(basiss,2)
#   R_t<-eval.penalty(basist,2)
#   
#   W_t_inv<-solve(W_t)
#   if(is.na(alpha_seq)[1]){
#     alpha_seq<-seq(-5,5)
#   }
#   
#   
#   inpr_vec<-list()
#   
#   parr_ridge<-function(ii){
#     
#     alpha<-alpha_seq[ii]
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     
#     for(ll in 1:K){
#       
#       Y_i<-Y_fd_cen[split_vec[[ll]]]
#       X_i<-X_new[split_vec[[ll]],]
#       X_minus<-X_new[-split_vec[[ll]],]
#       Y_minus<-Y_new[-split_vec[[ll]],]
#       
#       Beta_hat<-solve(t(X_minus)%*%X_minus+(10^alpha)*W_s)%*%t(X_minus)%*%Y_minus%*%W_t_inv
#       Y_hat<-fd(as.matrix(t(X_i%*%Beta_hat)),basist)
#       inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
#       
#     }
#     mean<-mean(unlist(inpr_vec))
#     sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
#     out<- list(mean=mean,
#                sd=sd)
#     return(out)
#   }
#   
#   
#   cores <- detectCores()
#   vec_par<-mclapply(seq(1,length(alpha_seq)),parr_ridge,mc.cores = cores)
#   par<-sapply(vec_par,"[[",1)
#   sd<-sapply(vec_par,"[[",2)
#   
#   l_opt<-alpha_seq[max(which(par==min(par)))]
#   
#   out<-list(par_opt=l_opt,
#             CV=par,
#             CV_sd=sd,
#             comb_list=alpha_seq,
#             X=X_fd,
#             Y=Y_fd,
#             type="RIDGE")
#   return(out)
# }
# 
# fr.ridge<-function(Y_fd,X_fd,basiss,basist,alpha_pen=10^0){
#   
#   X_fd_cen=center.fd(X_fd)
#   Y_fd_cen=center.fd(Y_fd)
#   X_mean<-mean.fd(X_fd)
#   Y_mean<-mean.fd(Y_fd)
#   n_obs<-length(X_fd$coefs[1,])
#   domain_s=X_fd$basis$rangeval
#   domain_t=Y_fd$basis$rangeval
#   
#   X_new<-inprod(X_fd_cen,basiss)
#   Y_new<-inprod(Y_fd_cen,basist)
#   
#   W_s<-eval.penalty(basiss)
#   W_t<-eval.penalty(basist)
#   
#   W_t_inv<-solve(W_t)
#   
#   cat("RIDGE:",alpha_pen,"     ")
#   B<-solve(t(X_new)%*%X_new+(alpha_pen)*W_s)%*%t(X_new)%*%Y_new%*%W_t_inv
#   Beta<-bifd(as.matrix(B),basiss,basist)
#   #Intercept
#   X_mean_new<-inprod(X_mean,basiss)
#   alpha<-Y_mean-fd(t(X_mean_new%*%B),basist)
#   
#   
#   out<-list(Beta_hat_fd=Beta,
#             B=B,
#             alpha=alpha,
#             alpha_pen=alpha_pen,
#             X=X_fd,
#             Y=Y_fd,
#             type="RIDGE")
#   
#   return(out)
# }
# 

# 
# 
# ## Historical ----------------------------------
# 
# 
# reg_his<-function(Y_fd,X_fd,K=10){
#   
#   n_basis_x_seq<-seq(20,40,by=5)
#   n_basis_y_seq<-seq(20,40,by=5)
#   n_obs<-length(X_fd$coefs[1,])
#   comb_list<-expand.grid(n_basis_x_seq,n_basis_y_seq)
#   inpr_vec<-numeric()
#   length_grid_int_t<-200
#   delta_t<-1/length_grid_int
#   grid_int_t<-seq(0,1,length.out = length_grid_int_t)
#   parr_his<-function(ii){
#     
#     basis<-as.numeric(comb_list[ii,])
#     n_basis_x<-basis[1]
#     n_basis_y<-basis[2]
#     basis_x <- create.bspline.basis(domain,nbasis = n_basis_x)
#     basis_y <- create.bspline.basis(domain,nbasis=n_basis_y)
#     
#     
#     in_prod_t<-function(t,X_fd){inprod(X_fd,basis_x,rng = c(0,t))}
#     X_t<-lapply(grid_int_t, in_prod_t,X_fd)
#     
#     
#     Y_eval<-eval.fd(Y_fd,grid_int_t)
#     basis_y_eval<-eval.basis(basis_y,grid_int_t)
#     eval_Y_t<-function(ii){(Y_eval[ii,]%o%basis_y_eval[ii,])}
#     Y_t<-lapply(1:length_grid_int_t,eval_Y_t)
#     eval_W_Y_t<-function(ii){(basis_y_eval[ii,]%o%basis_y_eval[ii,])}
#     W_Y_t<-lapply(1:length_grid_int_t,eval_W_Y_t )
#     
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     
#     for(ll in 1:K){
#       
#       X_i<-lapply(1:length_grid_int_t,function(ii){X_t[[ii]][split_vec[[ll]],]})
#       Y_fd_i<-Y_fd[split_vec[[ll]]]
#       X_minus<-lapply(1:length_grid_int_t,function(ii){X_t[[ii]][-split_vec[[ll]],]})
#       Y_minus<-lapply(1:length_grid_int_t,function(ii){Y_t[[ii]][-split_vec[[ll]],]})
#       
#       
#       A<-delta_t*Reduce("+",Map(function(ii){ kron(W_Y_t[[ii]],t(X_minus[[ii]])%*%X_minus[[ii]])},1:length_grid_int_t))
#       B<-vec(delta_t*Reduce("+",Map(function(ii){ t(X_minus[[ii]])%*%Y_minus[[ii]]},1:length_grid_int_t)))
#       B_his_vec<-ginv(A)%*%B
#       B_his<-matrix(B_his_vec,n_basis_x,n_basis_y)
#       aa<-lapply(1:length_grid_int_t,function(ii)X_i[[ii]]%*%B_his%*%basis_y_eval[ii,])
#       Y_hat_dis<-do.call(cbind,aa)
#       interp_basis<-create.bspline.basis(c(0,1),nbasis=length_grid_int_t,norder=1)
#       Y_hat<-fd(t(Y_hat_dis),interp_basis)
#       inpr_vec[ll]<-mean(diag(inprod(Y_fd_i-Y_hat,Y_fd_i-Y_hat)))
#     }
#     
#     (1/K)*sum(inpr_vec)
#     
#   }
#   
#   cores <- detectCores()
#   out<-mclapply(seq(1,length(comb_list[,1])),parr_his,mc.cores = cores)
#   
#   n_basis_opt<-as.numeric(comb_list[which(out==min(as.numeric(out))),])
#   
#   n_basis_x_opt<-n_basis_opt[1]
#   n_basis_y_opt<-n_basis_opt[2]
#   basis_x <- create.bspline.basis(domain,nbasis=n_basis_x_opt)
#   basis_y <- create.bspline.basis(domain,nbasis=n_basis_y_opt)
#   
#   in_prod_t<-function(t,X_fd){inprod(X_fd,basis_x,rng = c(0,t))}
#   X_t<-lapply(grid_int_t, in_prod_t,X_fd)
#   Y_eval<-eval.fd(Y_fd,grid_int_t)
#   basis_y_eval<-eval.basis(basis_y,grid_int_t)
#   eval_Y_t<-function(ii){(Y_eval[ii,]%o%basis_y_eval[ii,])}
#   Y_t<-lapply(1:length_grid_int_t,eval_Y_t)
#   eval_W_Y_t<-function(ii){(basis_y_eval[ii,]%o%basis_y_eval[ii,])}
#   W_Y_t<-lapply(1:length_grid_int_t,eval_W_Y_t )
#   
#   
#   A<-delta_t*Reduce("+",Map(function(ii){ kronecker(W_Y_t[[ii]],t(X_t[[ii]])%*%X_t[[ii]])},1:length_grid_int_t))
#   B<-vec(delta_t*Reduce("+",Map(function(ii){ t(X_t[[ii]])%*%Y_t[[ii]]},1:length_grid_int_t)))
#   B_his_vec<-ginv(A)%*%B
#   B_his<-matrix(B_his_vec,n_basis_x_opt,n_basis_y_opt)
#   
#   Beta_his<-bifd(B_his,basis_x,basis_y)
#   grid_int<-seq(0,1,length.out = 500)
#   
#   Beta_his_eval<-eval.bifd(grid_int,grid_int,Beta_his)
#   Beta_his_eval[lower.tri(Beta_his_eval,diag = T)]<-0
#   inter_basis<-create.bspline.basis(domain,nbasis = 500,norder = 1)
#   Beta_his<-bifd(Beta_his_eval,inter_basis,inter_basis)
#   
#   out<-list(B=B_his,
#             Beta_his=Beta_his,
#             basis_x_opt=basis_x,
#             basis_y_opt=basis_y,
#             n_basis_opt=n_basis_opt,
#             out=out)
#   
#   return(out)
#   
# }
# 
# ## Concurrent---------------------------------------
# 
# reg_conc<-function(Y_fd,X_fd,K=10){
#   
#   n_obs<-length(X_fd$coefs[1,])
#   lambda_vec<- 10^seq(-7,-4,by=1)
#   basis_con<-create.bspline.basis(c(0,1),nbasis = 70)
#   R<-eval.penalty(basis_con,2)
#   length_grid_int_t<-200
#   delta_t<-1/length_grid_int_t
#   grid_int_t<-seq(0,1,length.out = length_grid_int_t)
#   
#   
#   basis_con_eval<-t(eval.basis(grid_int_t,basis_con))
#   X_eval<-t(eval.fd(grid_int_t,X_fd))
#   in_prod_t<-function(ii){X_eval[,ii]%o%basis_con_eval[,ii]}
#   X_t<-lapply(1:length_grid_int_t, in_prod_t)
#   Y_eval<-t(eval.fd(grid_int_t,Y_fd))
#   
#   inpr_vec<-numeric()
#   parr_concu<-function(ii){
#     
#     lambda<-lambda_vec[ii]
#     ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
#     split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
#     
#     for(ll in 1:K){
#       
#       X_i<-lapply(1:length_grid_int_t,function(ii){X_t[[ii]][split_vec[[ll]],]})
#       Y_i<-Y_fd[split_vec[[ll]]]
#       X_minus<-lapply(1:length_grid_int_t,function(ii){X_t[[ii]][-split_vec[[ll]],]})
#       Y_minus<-Y_eval[-split_vec[[ll]],]
#       
#       A<-delta_t*Reduce("+",Map(function(ii){t(X_minus[[ii]])%*%X_minus[[ii]]},1:length_grid_int_t))
#       B<-delta_t*Reduce("+",Map(function(ii){t(X_minus[[ii]])%*%Y_minus[,ii]},1:length_grid_int_t))
#       B_his<-solve(A+lambda*R)%*%B
#       aa<-lapply(1:length_grid_int_t,function(ii)X_i[[ii]]%*%B_his)
#       Y_hat_dis<-do.call(cbind,aa)
#       interp_basis<-create.bspline.basis(c(0,1),nbasis=length_grid_int_t,norder=1)
#       Y_hat<-fd(t(Y_hat_dis),interp_basis)
#       
#       inpr_vec[ll]<-mean(tr(inprod(Y_i-Y_hat,Y_i-Y_hat)))
#     }
#     
#     (1/K)*sum(inpr_vec)
#     
#   }
#   
#   cores <- detectCores()
#   out<-mclapply(1:length(lambda_vec),parr_concu,mc.cores = cores)
#   l_opt<-lambda_vec[which(out==min(as.numeric(out)))]
#   
#   A<-delta_t*Reduce("+",Map(function(ii){t(X_t[[ii]])%*%X_t[[ii]]},1:length_grid_int_t))
#   B<-delta_t*Reduce("+",Map(function(ii){t(X_t[[ii]])%*%Y_eval[,ii]},1:length_grid_int_t))
#   B_his<-solve(A+l_opt*R)%*%B
#   
#   Beta_mon<-fd(B_his,basis_con)
#   grid_int<-seq(0,1,length.out = 500)
#   B_eval<-matrix(0,500,500)
#   diag(B_eval)<-eval.fd(grid_int,Beta_mon)
#   inter_basis<-create.bspline.basis(domain,nbasis =500,norder = 1)
#   Beta_con<-bifd(B_eval,inter_basis,inter_basis)
#   
#   out<-list(model=list(B=B,
#                        Beta_hat_fd=Beta_con),
#             lambda=l_opt)
#   return(out)
# }
# 
# ## ISE and PMSE ------------------------------------------------------------
# 
# get_ISE<-function(beta_hat_fd,Beta_vero_fd,case){
#   
#   length_grid_int<-500
#   delta<-1/length_grid_int
#   grid_int<-seq(0,1,length.out = length_grid_int)
#   
#   delta_bifd<-difference_bifd(Beta_vero_fd,beta_hat_fd)
#   eval_mat<-eval.bifd(grid_int,grid_int, delta_bifd)^2
#   
#   A<- eval.bifd(grid_int,grid_int,Beta_vero_fd)
#   
#   ind_0<-which(A==0,arr.ind = T)
#   sum_0<-sum(eval_mat[ind_0])
#   sum_1<-sum(eval_mat)-sum_0
#   area_0<-length(which(A==0))/(length_grid_int^2)
#   area_1<-1-area_0
#   ISE_0<-(1/area_0)*delta*delta*sum_0
#   ISE_1<-(1/area_1)*delta*delta*sum_1
#   
#   if(case=="Concurrent"){
#     ISE_1<-(sum(diag(eval_mat))*delta*delta)/(delta^2*length_grid_int)
#     ISE_tot<-sum(eval_mat)*delta*delta
#     ISE_0<-(ISE_tot-sum(diag(eval_mat))*delta*delta)/(1-delta^2*length_grid_int)
#   }
#   out<-list(ISE_0=ISE_0,
#             ISE_1=ISE_1)
#   return(out)
# }
# 
# difference_bifd<-function(bifd_1,bifd_2){
#   grid_s<-grid_t<-seq(0,1,length.out = 500)
#   X_1<-eval.bifd(grid_s,grid_t,bifd_1)
#   X_2<-eval.bifd(grid_s,grid_t,bifd_2)
#   diff<-(X_1-X_2)
#   B_spline_0<-create.bspline.basis(domain,nbasis = length(grid_s),norder = 1)
#   bifd(diff,B_spline_0,B_spline_0)
# }
# 
# get_PMSE<-function(Y_fd_test,X_fd_test,Beta,case_1=case){
#   
#   length_grid<-500
#   grid<-seq(0,1,length.out = length_grid)
#   delta<-1/length_grid
#   X_fd_eval<-t(eval.fd(grid,X_fd_test))
#   Y_fd_eval<-t(eval.fd(grid,Y_fd_test))
#   Beta_mat<-eval.bifd(grid,grid,Beta)
#   if(case_1=="Concurrent"){
#     Y_hat<-X_fd_eval%*%Beta_mat
#   }
#   else{
#     Y_hat<-delta*X_fd_eval%*%Beta_mat
#   }
#   PMSE<-mean(delta*rowSums((Y_fd_eval-Y_hat)^2))
#   
#   return(PMSE)
# }
# get_Area<-function(case,grid_int){
#   length_grid_int<-length(grid_int)
#   if(case=="Scenario I"){
#     Area_0<-1
#     Area_1<-0
#   }
#   if(case=="Scenario II"){
#     
#     Area_1<-0.25^2*pi
#     Area_0<-1-Area_1
#   }
#   if(case=="Scenario III"){
#     a=1
#     b=0.5
#     c=0.5
#     f_1<-function(s,t){((t-0.5)/c)^3+((t-0.5)/b)^2 + ((s-0.5)/a)^2-(0.5/a)^2}
#     f_2<-function(s,t){((t-0.5)/c)^3-((t-0.5)/b)^2 - ((s-0.5)/a)^2+(0.5/a)^2}
#     # tpos<-grid_int[((length_grid_int/2)+1):length_grid_int]
#     # tneg<-grid_int[1:(length_grid_int/2)]
#     z_pos<- outer(grid_int,grid_int,f_1)
#     ind_pos<-which(z_pos[,((length_grid_int/2)+1):length_grid_int]<0,arr.ind = T)
#     ind_pos[,2]<-ind_pos[,2]+(length_grid_int/2)
#     z_neg<- outer(grid_int,grid_int,f_2)
#     ind_neg<-which(z_neg[,1:(length_grid_int/2)]>0,arr.ind = T)
#     
#     ind_eval<-rbind(ind_neg,ind_pos)
#     
#     sum_int<-numeric(length(unique(ind_eval[,2])))
#     for (ii in 1:length(unique(ind_eval[,2]))) {
#       y<-unique(ind_eval[,2])[ii]
#       ind_1<-ind_eval[which(ind_eval[,2]==y),1]
#       
#       sum_int[ii]<-ind_1[length(ind_1)]-ind_1[1]
#     }
#     
#     delta=1/length_grid_int
#     Area_0<-delta^2*sum(sum_int)
#     Area_1<-1-Area_0
#   }
#   if(case=="Scenario IV"){
#     Area_0<-0
#     Area_1<-1
#   }
#   out<-list(Area_0=Area_0,
#             Area_1=Area_1)
#   return(out)
# }
# 
# mean_rep<-function (fdobj,nobs) 
# {
#   coef <- as.array(fdobj$coefs)
#   coefd <- dim(coef)
#   ndim <- length(coefd)
#   basisobj <- fdobj$basis
#   nbasis <- basisobj$nbasis
#   coefmean <- apply(coef, 1, mean)
#   mean_rep <- fd(matrix(rep(coefmean,nobs),coefd[1],nobs), basisobj)
#   return(mean_rep)
# }
# innerProd<-function (A, B, gridbrk, intlen) 
# {
#   brkX <- gridbrk
#   oupmat <- array(0, c(dim(A)[1], dim(B)[2]))
#   for (i in 1:(length(brkX) - 1)) {
#     a1 <- brkX[i] + 1
#     a2 <- brkX[i + 1]
#     A[, a1] <- A[, a1]/2
#     A[, a2] = A[, a2]/2
#     oupmat <- oupmat + A[, a1:a2] %*% B[a1:a2, ] * intlen[i]/(a2 - 
#                                                                 a1)
#   }
#   return(oupmat)
# }
# ## TNR and FNR-----
# get_TNR_FNR<-function(beta_hat_fd,Beta_vero_fd){
#   
#   length_grid_int<-500
#   delta<-1/length_grid_int
#   grid_int<-seq(0,1,length.out = length_grid_int)
#   
#   eval_beta_vero<-eval.bifd(grid_int,grid_int, Beta_vero_fd)
#   eval_beta_hat<-eval.bifd(grid_int,grid_int, beta_hat_fd)
#   area_0<-length(which(eval_beta_vero==0))/(length_grid_int^2)
#   area_1<-1-area_0
#   ind_0<-which(eval_beta_vero==0)
#   TNR<-length(which(eval_beta_hat[ind_0]==0))/length(ind_0)
#   ind_1<-which(eval_beta_vero!=0)
#   FNR<-length(which(eval_beta_hat[ind_1]==0))/length(ind_1)
#   
#   out<-list(TNR=TNR,
#             FNR=FNR)
#   return(out)
# }
# 
# R2<-function(mod){
#   Y<-mod$Y
#   X<-mod$X
#   Y=center.fd(Y)
#   X=center.fd(X)
#   rangevals<-mod$Beta_hat_fd$sbasis$rangeval
#   rangevalt<-mod$Beta_hat_fd$tbasis$rangeval
#   X_new=inprod(X,mod$Beta_hat_fd$sbasis)
#   Y_hat=fd(t(X_new%*%mod$B),mod$Beta_hat_fd$tbasis)
#   
#   
#   
#   length_grid<-500
#   
#   grid_s<-seq(rangevals[1],rangevals[2],length.out = length_grid)
#   grid_t<-seq(rangevalt[1],rangevalt[2],length.out = length_grid)
#   
#   X_mat<-t(eval.fd(grid_s,X))
#   Y_mat<-t(eval.fd(grid_t,Y))
#   Y_hat_mat<-t(eval.fd(grid_t,Y_hat))
#   
#   SSres<-colSums((Y_mat-Y_hat_mat)^2)
#   SStot<-colSums((Y_mat)^2)
#   ss_div<-SSres/SStot
#   ss_div[which(ss_div<=0)]<-0
#   ss_div[which(ss_div>=1)]<-1
#   R2t<-fd(1-ss_div,create.bspline.basis(c(rangevalt[1],rangevalt[2]),length_grid))
#   
#   R2<-1/(rangevalt[2]-rangevalt[1])* inprod(R2t,fd(1,create.constant.basis(rangevalt)))
#   
#   out<-list(R2=R2,
#             R2t=R2t)
#   return(out)
# }
# get_per_0<-function(mod){
#   length_grid<-500
#   if(is.null(mod$Beta_hat_fd$sbasis$rangeval)){
#     rangevals<-mod$sbasis$rangeval
#     rangevalt<-mod$tbasis$rangeval
#     beta<-mod
#   }
#   else{
#     rangevals<-mod$Beta_hat_fd$sbasis$rangeval
#     rangevalt<-mod$Beta_hat_fd$tbasis$rangeval
#     beta<-mod$Beta_hat_fd
#   }
#   grid_s<-seq(rangevals[1],rangevals[2],length.out = length_grid)
#   grid_t<-seq(rangevalt[1],rangevalt[2],length.out = length_grid)
#   A=eval.bifd(seq(rangevals[1],rangevals[2],length.out = 500),seq(rangevalt[1],rangevalt[2],length.out = 500),beta)
#   out<-length(which(A==0))/(500*500)
#   return(out)
# }
# # 
# #   matplot(t((1/1000)*t(eval.fd(grid_int,X_fd_tra[split_vec[[ll]]]))%*%eval.bifd(grid_int,grid_int,Beta_vero_fd)),type = "l")
# #   matplot(t((1/1000)*t(eval.fd(grid_int,X_fd_tra[split_vec[[ll]]]))%*%eval.bifd(grid_int,grid_int,mod_ridge$model$Beta_fd)),type = "l")
# #   matplot(t((1/1000)*t(eval.fd(grid_int,X_fd_tra[split_vec[[ll]]]))%*%eval.bifd(grid_int,grid_int,Beta_cpp)),type = "l")
# #   x11()
# #   par(mfrow=c(1,3))
# #   plot(eval.bifd(grid_int,0.5,Beta_vero_fd),ylim = c(0,1))
# #   plot(eval.bifd(grid_int,0.5,Beta_cpp),ylim = c(0,1))
# #   plot(eval.bifd(grid_int,0.5,mod_ridge$model$Beta_fd),ylim = c(0,1))
# #   delta*(sum((eval.bifd(grid_int[1:100],0.5,mod_ridge$model$Beta_fd))^2)+sum((eval.bifd(grid_int[401:500],0.5,mod_ridge$model$Beta_fd)^2)))
# #   delta*(sum((eval.bifd(grid_int[1:100],0.5,Beta_cpp))^2)+sum((eval.bifd(grid_int[401:500],0.5,Beta_cpp)^2)))
# #   
# #   
# #   (sum(eval_mat[1:100,250])+sum(eval_mat[401:500,250]))*delta
# #   (sum(eval_mat[101:400,250]))*delta
# 
# 
