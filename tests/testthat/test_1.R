


case<-"Scenario III"
n_obs_tra<-150

n_obs<-n_obs_tra

data<-simulate_data(case,n_obs=n_obs)
Beta_vero_mat<-data$beta_matrix
X_fd=data$X_fd
Y_fd=data$Y_fd
# Basis s and t
domain=c(0,1)
n_basis_s<-60
n_basis_t<-60
breaks_s<-seq(0,1,length.out = (n_basis_s-2))
breaks_t<-seq(0,1,length.out = (n_basis_t-2))
basis_s <- fda::create.bspline.basis(domain,breaks=breaks_s)
basis_t <- fda::create.bspline.basis(domain,breaks=breaks_t)


# 
 mod_slasso_cv<-slasso.fr_cv(Y_fd = Y_fd,X_fd=X_fd,basis_s=basis_s,basis_t=basis_t,lambdas_L=seq(0,1,by=1),lambdas_s=c(-9),lambdas_t=-7,B0=NULL,lam_opt_method="min",
                            ss_rule_par=c(0.1,0.5,1),max_iterations=1000,X_test = X_fd_test[1:500],Y_test = Y_fd_test[1:500],K=2,invisible=1,cores=2)
mod_slasso<-slasso.fr(Y_fd = Y_fd,X_fd=X_fd,basis_s=basis_s,basis_t=basis_t,lambdas_L = -1.5,lambdas_s =-8,lambdas_t = -7,B0 =NULL,invisible=1,max_iterations=10)

slasso::plot.slasso_cv(mod_slasso_cv)
slasso::plot.slasso(mod_slasso)
