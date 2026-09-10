# Data-dependent rebuilding only. No method packages are loaded.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
input<-normalizePath(args[1]);out<-args[2];stopifnot(!file.exists(out));dir.create(out,recursive=TRUE)
source('study/lib/reference_endpoints.R')
math<-re_load_math('study/lib/endpoint_math.R',re_sha('study/lib/endpoint_math.R'))
units<-re_read(file.path(input,'units.tsv'));stopifnot(nrow(units)==120L)
reads<-list();options(digits=17)
for(i in seq_len(nrow(units))){
 p<-file.path(input,paste0(units$task_id[i],'.rds'));z<-readRDS(p)
 stopifnot(z$schema=='PUBLIC_SIMULATION_EVIDENCE_V1',z$unit$task_id==units$task_id[i])
 result<-re_unit_endpoints(z$unit,z$genes,z$cells,z$samples,z$evaluation,z$regimes,math)
 for(n in names(result)){f<-file.path(out,paste0(n,'.tsv'));re_write(result[[n]],f,file.exists(f))}
 reads[[i]]<-data.frame(path=paste0('simulation/',basename(p)),role='PUBLIC_FROZEN_NUMERIC_BUNDLE',bytes=file.info(p)$size,sha256=re_sha(p))
 cat('REFERENCE_REBUILD',i,'/120\n');flush.console()
}
re_write(do.call(rbind,reads),file.path(out,'ACTUAL_ENDPOINT_READS.tsv'))
writeLines('DERIVED_ONLY_NO_MODELS_NO_NNLS',file.path(out,'DERIVED_BUILD_COMPLETE'))
