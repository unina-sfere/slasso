
data<-simulate_data("Scenario II",n_obs=1000)
X_fd=data$X_fd
Y_fd=data$Y_fd
domain=c(0,1)
n_basis_s<-30
n_basis_t<-30
breaks_s<-seq(0,1,length.out = (n_basis_s-2))
breaks_t<-seq(0,1,length.out = (n_basis_t-2))
basis_s <- fda::create.bspline.basis(domain,breaks=breaks_s)
basis_t <- fda::create.bspline.basis(domain,breaks=breaks_t)


mod_smooth_cv<-fr.usc.cv(Y_fd_tra,X_fd_tra,basis_x,basis_y,K=10,lambdas_s = 10^seq(-6,-2),lambdas_t = 10^seq(-6,-2))
mod_smooth<-fr.usc(Y_fd,X_fd,basis_s,basis_t,lambdas_s=10^-6,lambdas_t =10^-6)
lambda_L_vec=10^seq(0,1,by=0.1) 
lambda_s_vec=10^seq(-6,-5) 
lambda_t_vec=10^seq(-5,-5) 
mod_slasso_cv<-slasso.fr_cv(Y_fd = Y_fd,X_fd=X_fd,basis_s=basis_s,basis_t=basis_t,lambda_L_vec = lambda_L_vec,lambda_s_vec = lambda_s_vec,lambda_t_vec =lambda_t_vec,max_iterations=1000,K=10,invisible=1,ncores=12)
mod_slasso<-slasso.fr(Y_fd = Y_fd,X_fd=X_fd,basis_s=basis_s,basis_t=basis_t,lambda_L = mod_slasso_cv$lambda_opt_vec[1],lambda_s = mod_slasso_cv$lambda_opt_vec[2],lambda_t =  mod_slasso_cv$lambda_opt_vec[3],invisible=1,max_iterations=1000)
mod_slasso<-slasso.fr(Y_fd = Y_fd,X_fd=X_fd,basis_s=basis_s,basis_t=basis_t,lambda_L = 3.16,lambda_s = 10^-5,lambda_t =  10^-5,B0 =NULL,invisible=1,max_iterations=10000)

 plot(mod_slasso_cv)
plot(mod_slasso)
