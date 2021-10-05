# C++ objective and gradient functions---------------------------------------------------------------------

objective.include <- '
Rcpp::NumericVector objective_fun(SEXP B_mat,SEXP env) {

/*B=matrix(B,nrow=n_basis_x,ncol=n_basis_y,byrow = FALSE)
sum(diag( t(Y_coef_tra-X_coef_tra%*%B)%*%(Y_coef_tra-X_coef_tra%*%B)))+lambda*delta_x^2*sum(abs(t(basis_values_x)%*%B%*%basis_values_y))*/

Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
int n_basis_x = Rcpp::as<int>(e["n_basis_x"]);
int n_basis_y = Rcpp::as<int>(e["n_basis_y"]);
double lambda_x_opt=Rcpp::as<double>(e["lambda_x_opt"]);
double lambda_y_opt=Rcpp::as<double>(e["lambda_y_opt"]);


Rcpp::NumericVector B_nv =Rcpp::as<Rcpp::NumericVector>(B_mat);

arma::mat B_vec=Rcpp::as<Rcpp::NumericVector>(B_nv);

arma::mat Y_newc= Rcpp::as<arma::mat>(e["Y_newc"]);
arma::mat X_newc = Rcpp::as<arma::mat>(e["X_newc"]);


arma::mat W_X = Rcpp::as<arma::mat>(e["W_X"]);
arma::mat W_Y = Rcpp::as<arma::mat>(e["W_Y"]);
arma::mat R_X = Rcpp::as<arma::mat>(e["R_X"]);
arma::mat R_Y = Rcpp::as<arma::mat>(e["R_Y"]);

arma::mat B(n_basis_x,n_basis_y);

B=arma::reshape(B_vec,n_basis_x,n_basis_y);

Rcpp::NumericVector out =Rcpp::as<Rcpp::NumericVector>(wrap(trace(-2*trans(Y_newc)*X_newc*B)+trace(W_Y*trans(B)*trans(X_newc)*X_newc*B)+trace(W_Y*trans(Y_newc)*Y_newc)
+lambda_x_opt*trace(trans(B)*R_X*B*W_Y)+lambda_y_opt*trace(trans(B)*W_X*B*R_Y)));

return out;
}'



gradient.include <- 'Rcpp::NumericVector grad_fun(SEXP B_mat, SEXP env) {
/* B<-matrix(B,nrow=n_basis_x,ncol=n_basis_y,byrow = FALSE)
f<-function(x,y){as.numeric(sign(t(x)%*%B%*%y))*x%o%y}
sum<-Reduce("+",Map(f,basis_values_x_list,basis_values_y_list))
c(-2*t(X_coef_tra)%*%Y_coef_tra+2*t(X_coef_tra)%*%X_coef_tra%*%B+lambda* delta_x^2*sum) */


Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
int n_basis_x = Rcpp::as<int>(e["n_basis_x"]);
int n_basis_y = Rcpp::as<int>(e["n_basis_y"]);

double lambda_x_opt=Rcpp::as<double>(e["lambda_x_opt"]);
double lambda_y_opt=Rcpp::as<double>(e["lambda_y_opt"]);


Rcpp::NumericVector B_nv =Rcpp::as<Rcpp::NumericVector>(B_mat);

arma::mat B_vec=Rcpp::as<Rcpp::NumericVector>(B_nv);

arma::mat Y_newc = Rcpp::as<arma::mat>(e["Y_newc"]);
arma::mat X_newc = Rcpp::as<arma::mat>(e["X_newc"]);

arma::mat  W_X = Rcpp::as<arma::mat >(e["W_X"]);
arma::mat  W_Y = Rcpp::as<arma::mat >(e["W_Y"]);
arma::mat  R_X = Rcpp::as<arma::mat >(e["R_X"]);
arma::mat  R_Y = Rcpp::as<arma::mat >(e["R_Y"]);

arma::mat B(n_basis_x,n_basis_y);

B=arma::reshape(B_vec,n_basis_x,n_basis_y);

Rcpp::NumericVector out =Rcpp::as<Rcpp::NumericVector>(wrap(-2*trans(X_newc)*Y_newc+2*trans(X_newc)*X_newc*B*W_Y
 +lambda_x_opt*2*R_X*B*W_Y+lambda_y_opt*2*W_X*B*R_Y));

return out;
}'


objective.body <- '
typedef Rcpp::NumericVector (*funcPtr)(SEXP,SEXP);
return(XPtr<funcPtr>(new funcPtr(&objective_fun)));'
gradient.body <- '
typedef Rcpp::NumericVector (*funcPtr)(SEXP,SEXP);
return(XPtr<funcPtr>(new funcPtr(&grad_fun)));'

objective <- cxxfunction(signature(), body=objective.body,
                         inc=objective.include, plugin="RcppArmadillo")

gradient <- cxxfunction(signature(), body=gradient.body,
                        inc=gradient.include, plugin="RcppArmadillo")




# Kronecker ---------------------------------------------------------------

kron <- cxxfunction(signature(A="SEXP",B="SEXP"), body='
                    arma::mat A_m = Rcpp::as<arma::mat>(A);
                    arma::mat B_m = Rcpp::as<arma::mat>(B);
                    arma::mat z = kron(A_m,B_m);
                    return wrap( z );',
                    plugin="RcppArmadillo")


# Generalized Inverse -----------------------------------------------------------------

ginverse <- cxxfunction(signature(A="SEXP"), body='
                        #define ARMA_64BIT_WORD 1
                        #include <RcppArmadillo.h>
                        // [[Rcpp::depends(RcppArmadillo)]]
                        // [[Rcpp::plugins(cpp11)]] 
                        arma::mat A_m = Rcpp::as<arma::mat>(A);
                        arma::mat z = pinv(A_m);
                        return wrap( z );',
                        plugin="RcppArmadillo")




# P function for approximated df ---------------------------------------------------------------

P_func <- cxxfunction(signature(A="SEXP",B="SEXP",B_l="SEXP"), body='
                      arma::mat basis_values_x = Rcpp::as<arma::mat>(A);
                      arma::mat basis_values_y = Rcpp::as<arma::mat>(B);
                      arma::mat P_m = Rcpp::as<arma::mat>(B_l);
                      arma::mat somma(basis_values_x.n_rows*basis_values_y.n_rows,basis_values_x.n_rows*basis_values_y.n_rows);
                      somma.zeros(basis_values_x.n_rows*basis_values_y.n_rows,basis_values_x.n_rows*basis_values_y.n_rows);
                      double lenght_grid=basis_values_x.n_cols;
                      arma::mat D(lenght_grid,lenght_grid);
                      D.zeros(lenght_grid,lenght_grid);
                      for (int i=0; i<lenght_grid ;i++){
                      D.diag()=P_m.row(i);
                      somma= somma + kron(basis_values_y*D*trans(basis_values_y),basis_values_x.col(i)*trans(basis_values_x.col(i)));
                      Rcpp::Rcout << i;
                      }
                      return wrap( somma );',
                      plugin="RcppArmadillo")

