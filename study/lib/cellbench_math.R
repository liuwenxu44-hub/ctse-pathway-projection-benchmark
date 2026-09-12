cor_reason <- function(x,y) {
 stopifnot(length(x)==length(y));ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
 if(length(x)<3)return('CORRELATION_UNDEFINED_INSUFFICIENT_PAIRS')
 a<-sd(x)<=1e-12*max(1,max(abs(x)));b<-sd(y)<=1e-12*max(1,max(abs(y)))
 if(a&&b)return('CORRELATION_UNDEFINED_ZERO_VARIANCE_BOTH')
 if(a)return('CORRELATION_UNDEFINED_ZERO_VARIANCE_ESTIMATE')
 if(b)return('CORRELATION_UNDEFINED_ZERO_VARIANCE_REFERENCE')
 'CORRELATION_DEFINED'
}
safe_cor <- function(x,y,method='spearman') {
 if(cor_reason(x,y)!='CORRELATION_DEFINED')return(NA_real_)
 ok<-is.finite(x)&is.finite(y);unname(cor(x[ok],y[ok],method=method))
}
fraction <- function(x,y) {
 if(length(x)!=3||length(y)!=3||any(!is.finite(c(x,y))))return(c(MAE=NA_real_,RMSE=NA_real_))
 c(MAE=mean(abs(x-y)),RMSE=sqrt(mean((x-y)^2)))
}
sne <- function(x,y) if(length(x)!=length(y)||length(y)<2||any(!is.finite(c(x,y)))) NA_real_ else sqrt(mean((x-y)^2))/max(sd(y),1e-8)
projection <- function(delta,idx) {
 if(any(!is.finite(delta))||!length(idx))return(NA_real_)
 z<-sqrt(sum(delta^2));if(z<=1e-15)return(0)
 sum(delta[idx])/(sqrt(length(idx))*z)
}
amplitude <- function(delta,anchor) {
 if(any(!is.finite(c(delta,anchor))))return(NA_real_)
 z<-sqrt(sum(anchor^2));if(z<=0)return(NA_real_)
 sqrt(sum(delta^2))/z
}
two_mean <- function(x,strata) {z<-vapply(split(x,strata),function(y)if(any(is.finite(y)))mean(y[is.finite(y)])else NA_real_,0.0);if(any(is.finite(z)))mean(z[is.finite(z)])else NA_real_}
vecmean <- function(x,strata) {
 z<-lapply(split(seq_len(ncol(x)),strata),function(i) {q<-x[,i,drop=FALSE];q<-q[,colSums(!is.finite(q))==0,drop=FALSE];if(!ncol(q))return(NULL);rowMeans(q)})
 z<-Filter(Negate(is.null),z);if(!length(z))return(rep(NA_real_,nrow(x)))
 rowMeans(do.call(cbind,z))
}
direction <- function(x) ifelse(!is.finite(x),NA_character_,ifelse(x< -0.05,'negative',ifelse(x>0.05,'positive','zero')))
write_tsv <- function(x,p) {old<-options(digits=17);on.exit(options(old));write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')}
