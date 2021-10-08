#' @title Plot the results of  the S-LASSO method
#' @description This function provides plots of the estimated cluster mean functions and of the classified curves when applied to the output of `sasfclust`, whereas
#'  provides the cross-validation plots when applied to the output of `sasfclust_cv`. In the latter case the first plot displays the CV values as a function of  \code{G}, \code{lambda_s} and \code{lambda_l};
#'  the second plot displays the CV values as a function of \code{lambda_s} and \code{lambda_l} for \code{G} fixed at its optimal value;
#'   the third plot displays the CV values as a function of \code{lambda_l} for \code{G} and \code{lambda_s}   fixed at their optimal value.
#'
#' @param x The output of  either `slasso.fr_cv` or `slasso.fr`.
#' @param ... No additional parameters, called for side effects.
#' @return No return value, called for side effects.
#' @rdname plot.slasso
#' @method plot slasso_cv
#' @export
#' @examples
#' \donttest{
#' library(slasso)
#' 
#' }
plot.slasso_cv<-function(mod,type="",ss_rule_par=1,fun_c="mean",ylim=NA,...){
  
  if(fun_c=="min"){
    fun=colMins
  }
  else if(fun_c=="mean"){
    fun=colMeans
  }
  
  
  if(mod$type!="SLASSO_CV") stop("Wrong model provided!")
  x<-seq(1,length(mod$comb_list[,1]))
  CV_i=mod$CV
  sd_i=mod$CV_sd
  zeros_i=mod$per_0
  labels_L<-lapply(1:length(mod$comb_list[,1]),function(ii){a<-as.character(signif(mod$comb_list[ii,],digits = 1));paste(a[1],a[2],a[3])})
  
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit( graphics::par(oldpar))
  graphics::layout(matrix(rbind(c(1,1,1),c(3,2,3)),2,3))
  graphics::par(mar=c(6,6,6,5))
  base::plot(CV_i,pch=16,cex=0.5,col=2,type="l",xaxt="n",xlab="",ylim=c(min(CV_i)-1.01*max(sd_i),max(CV_i)+1.01*max(sd_i)),ylab="CV")
  graphics::points(CV_i,pch=16,cex=0.5,col=2)
  graphics::segments(x-0.1,CV_i+sd_i,x+0.1)
  graphics::segments(x-0.1,CV_i-sd_i,x+0.1)
  graphics::mtext(text=labels_L,side=1,at=x,las=2,cex=0.75)
  graphics::mtext(text=as.character(round(zeros_i*100)),side=3,at=x,las=2,cex=0.75)
  # abline(v=which(CV_i==min(CV_i)))
  abline(h=min(CV_i))
  lamb_s<-unique(mod$comb_list[,2])
  lamb_t<-unique(mod$comb_list[,3])
  lamb_L<-unique(mod$comb_list[,1])
  len_s<-length(lamb_s)
  len_t<-length(lamb_t)
  len_L<-length(lamb_L)
  
  
  
  ###Lambda L
  mat_CV<-matrix(0,len_s*len_t,len_L)
  mat_CV_sd<-matrix(0,len_s*len_t,len_L)
  mat_per_0<-matrix(0,len_s*len_t,len_L)
  for(ii in 1:len_L){
    mat_CV[,ii]<-CV_i[which(mod$comb_list[,1]==lamb_L[ii])]
    mat_CV_sd[,ii]<-sd_i[which(mod$comb_list[,1]==lamb_L[ii])]
    mat_per_0[,ii]<-mod$per_0[which(mod$comb_list[,1]==lamb_L[ii])]
  }
  
  min_CV_L<-fun(mat_CV)
  aa<-sapply(1:len_L, function(ii)which(CV_i==min_CV_L[ii]))
  sd_CV_L<-sd_i[aa]
  min_per_0<-zeros_i[aa]
  labels_L<-labels_L[aa]
  
  x<-log10(lamb_L)
  base::plot(x,min_CV_L,pch=16,cex=0.5,col=2,type="l",xaxt="n",xlab="",ylim=c(min(min_CV_L)-1.01*max(sd_CV_L),max(min_CV_L)+1.01*max(sd_CV_L)),ylab="CV in function of lambda_L")
  graphics::points(min_CV_L,pch=16,cex=0.5,col=2)
  graphics::segments(x-0.1,min_CV_L+sd_CV_L,x+0.1)
  graphics::segments(x-0.1,min_CV_L-sd_CV_L,x+0.1)
  graphics::mtext(text=labels_L,side=1,at=x,las=2,cex=0.75)
  graphics::mtext(text=as.character(round(min_per_0*100)),side=3,at=x,las=2,cex=0.75)
  graphics::abline(h=min(min_CV_L),col=2,lty=2)
  
  
  
}


#' @rdname plot.slasso
#' @method plot slasso
#' @export
#'
plot.slasso<-function(mod,...){
  
  length_grid=200
  rangevals<-mod$Beta_hat_fd$sbasis$rangeval
  rangevalt<-mod$Beta_hat_fd$tbasis$rangeval
  A=fda::eval.bifd(seq(rangevals[1],rangevals[2],length.out = length_grid),seq(rangevalt[1],rangevalt[2],length.out = length_grid),mod$Beta_hat_fd)
  names(A)<-c("s","t")
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit( graphics::par(oldpar))
  graphics::par(mfrow=c(1,2))
  plot3D::contour2D(z=A,x=seq(rangevals[1],rangevals[2],length.out = length_grid),y=seq(rangevalt[1],rangevalt[2],length.out = length_grid),xlab="s",ylab="t")
  persp(A,x=seq(rangevals[1],rangevals[2],length.out = length_grid),y=seq(rangevalt[1],rangevalt[2],length.out = length_grid),zlab="",xlab="s",ylab="t",ticktype="detailed",col="lightblue",shade = 0.75,border = NA)
  
}
