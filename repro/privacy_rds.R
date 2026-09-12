# Inspect every character field and attribute, without modifying the objects.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
root<-normalizePath(args[1]);out<-args[2];stopifnot(!file.exists(out))
patterns<-c('/(home|Users|mnt)/','/data/[A-Za-z][A-Za-z0-9_-]*/',
            '[A-Z]:[\\\\/](Users|work|proj)',
            '172[.](1[6-9]|2[0-9]|3[01])[.][0-9]+[.][0-9]+',
            'github_pat_[A-Za-z0-9_]{20,}','gh[pousr]_[A-Za-z0-9]{20,}')
check<-function(x){
 if(is.environment(x)||is.function(x)||typeof(x)%in%c('externalptr','weakref'))return(FALSE)
 if(is.character(x)&&any(vapply(patterns,function(p)any(grepl(p,x,perl=TRUE),na.rm=TRUE),logical(1))))return(FALSE)
 if(is.list(x)&&!all(vapply(x,check,logical(1))))return(FALSE)
 att<-attributes(x)
 if(!is.null(att)&&!all(vapply(att,check,logical(1))))return(FALSE)
 TRUE
}
files<-list.files(root,pattern='[.]rds$',recursive=TRUE,full.names=TRUE);rows<-list()
for(i in seq_along(files)){
 p<-files[i];ok<-check(readRDS(p))
 rows[[i]]<-data.frame(file=substring(p,nchar(root)+2L),status=if(ok)'PASS'else'PRIVATE_OR_EXECUTABLE_CONTENT')
 if(i%%100L==0L){cat('RDS_PRIVACY_CHECKED',i,'/',length(files),'\n');flush.console()}
}
write.table(do.call(rbind,rows),out,sep='\t',quote=FALSE,row.names=FALSE)
stopifnot(all(vapply(rows,function(z)z$status=='PASS',logical(1))))
cat('PUBLIC_RDS_PRIVACY_PASS',length(files),'\n')
