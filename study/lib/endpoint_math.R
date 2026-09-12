safe_mean<-function(x)if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
safe_rank<-function(x,y)if(length(x)>=3L&&sd(x)>0&&sd(y)>0)cor(x,y,method='spearman') else NA_real_
classify<-function(x)ifelse(!is.finite(x),NA_character_,ifelse(abs(x)<=.05,'zero',ifelse(x>0,'positive','negative')))
projection<-function(delta,genes,pathway){
 set<-genes%in%sprintf('sg%04d',if(pathway=='signal')1:50 else 51:100)
 norm<-sqrt(sum(delta^2));n<-sum(set)
 if(!n)return(c(P=NA_real_,norm=norm,members=0))
 c(P=if(norm<=1e-15)0 else sum(delta[set])/(sqrt(n)*norm),norm=norm,members=n)
}
strata_masks<-function(ev,genes,k){
 changed<-ev$design_changed[genes,k]
 list(ALL=rep(TRUE,length(genes)),DIRECT_SIGNAL=genes%in%sprintf('sg%04d',1:50),NEGATIVE_CONTROL=genes%in%sprintf('sg%04d',51:100),BACKGROUND=genes%in%sprintf('sg%04d',501:2000),DESIGN_CHANGED=changed,DESIGN_UNCHANGED=!changed)
}
class_summary<-function(truth,prediction){
 cl<-c('negative','zero','positive');ok<-!is.na(prediction)
 rec<-vapply(cl,function(z){j<-truth==z&ok;if(any(j))mean(prediction[j]==z) else NA_real_},numeric(1))
 nz<-truth!='zero';zero<-!nz
 c(planned=length(truth),evaluable=sum(ok),coverage=mean(ok),overall_accuracy=if(any(ok))mean(truth[ok]==prediction[ok]) else NA_real_,balanced_accuracy_macro_recall=if(all(is.finite(rec)))mean(rec) else NA_real_,nonzero_planned=sum(nz),nonzero_evaluable=sum(nz&ok),nonzero_sign_accuracy=if(any(nz&ok))mean(prediction[nz&ok]==truth[nz&ok]) else NA_real_,zero_planned=sum(zero),zero_evaluable=sum(zero&ok),zero_truth_recall=if(any(zero&ok))mean(prediction[zero&ok]=='zero') else NA_real_,recall_negative=rec[1],recall_positive=rec[3])
}
