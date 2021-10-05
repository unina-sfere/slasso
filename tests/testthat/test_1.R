


case<-"Scenario III"
n_obs_tra<-150
n_obs_test<-4000
n_obs<-n_obs_tra+n_obs_test

data<-simulate_data(case,n_obs=n_obs)
Beta_vero_mat<-data$beta_matrix
X<-data$X
Y<-data$Y

domain<-c(0,1)
length_grid<-500
norder<-4
grid<-seq(0,1,length.out = length_grid)


# Smoothing
n_basis_x<-min(80,length_grid)
n_basis_y<-min(80,length_grid)
breaks_x<-seq(0,1,length.out = (n_basis_x-2))
breaks_y<-seq(0,1,length.out = (n_basis_y-2))
basis_x <- create.bspline.basis(domain,breaks=breaks_x)
basis_y <- create.bspline.basis(domain,breaks=breaks_y)
X_fd <- smooth.basis(grid,X,basis_x)$fd
Y_fd <- smooth.basis(grid,Y,basis_y)$fd
inter_basis<-create.bspline.basis(domain,nbasis = length(grid),norder = 1)
Beta_vero_fd<-bifd(Beta_vero_mat,inter_basis,inter_basis)

# Basis s and t
n_basis_x<-min(60,length_grid)
n_basis_y<-min(60,length_grid)
breaks_x<-seq(0,1,length.out = (n_basis_x-2))
breaks_y<-seq(0,1,length.out = (n_basis_y-2))
basis_x <- create.bspline.basis(domain,breaks=breaks_x)
basis_y <- create.bspline.basis(domain,breaks=breaks_y)

# Matrices
W_X<-eval.penalty(basis_x)
W_X_sp<-sparseMatrix(which(W_X>0,arr.ind=T)[,1],which(W_X>0,arr.ind=T)[,2],x=W_X[which(W_X>0,arr.ind=T)])
W_Y<-eval.penalty(basis_y)
W_Y_sp<-sparseMatrix(which(W_Y>0,arr.ind=T)[,1],which(W_Y>0,arr.ind=T)[,2],x=W_Y[which(W_Y>0,arr.ind=T)])
W_XY<-inprod(basis_x,basis_y)
W_XY_sp<-sparseMatrix(which(W_XY>0,arr.ind=T)[,1],which(W_XY>0,arr.ind=T)[,2],x=W_XY[which(W_XY>0,arr.ind=T)])
R_X<-eval.penalty(basis_x,2)
R_X_sp<-sparseMatrix(which(R_X!=0,arr.ind=T)[,1],which(R_X!=0,arr.ind=T)[,2],x=R_X[which(R_X!=0,arr.ind=T)])
R_Y<-eval.penalty(basis_y,2)
R_Y_sp<-sparseMatrix(which(R_Y!=0,arr.ind=T)[,1],which(R_Y!=0,arr.ind=T)[,2],x=R_Y[which(R_Y!=0,arr.ind=T)])

# Training and Test set
X_fd_tra_nc<-X_fd[1:n_obs_tra]
Y_fd_tra_nc<-Y_fd[1:n_obs_tra]
X_fd_tra<-center.fd(X_fd_tra_nc)
Y_fd_tra<-center.fd(Y_fd_tra_nc)
X_fd_test<-X_fd[(n_obs_tra+1):(n_obs)]-mean_rep(X_fd_tra_nc,n_obs_test)
Y_fd_test<-Y_fd[(n_obs_tra+1):(n_obs)]-mean_rep(Y_fd_tra_nc,n_obs_test)
X_coef_tra<-t(X_fd_tra$coefs)
Y_coef_tra<-t(Y_fd_tra$coefs)
X_coef_test<-t(X_fd_test$coefs)
Y_coef_test<-t(Y_fd_test$coefs)   



# 
 mod_slasso_cv<-slasso.fr.cv(Y_fd = Y_fd_tra,X_fd=X_fd_tra,basis_s=basis_x,basis_t=basis_y,lambdas_L=seq(0,1,by=1),lambdas_s=c(-9),lambdas_t=-7,B0=NULL,lam_opt_method="min",
                            ss_rule_par=c(0.1,0.5,1),max_iterations=1000,X_test = X_fd_test[1:500],Y_test = Y_fd_test[1:500],K=2,invisible=1,cores=2)
mod_slasso<-slasso.fr(Y_fd_tra,X_fd_tra,basis_x,basis_y,lambdas_L = -1.5,lambdas_s =-8,lambdas_t = -7,B0 =NULL,invisible=1,max_iterations=10)
