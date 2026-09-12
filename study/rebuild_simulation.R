# Original positive-control support, separate from paired-reference sensitivity.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
input<-normalizePath(args[1]);out<-args[2];stopifnot(!file.exists(out));dir.create(out,recursive=TRUE)
source('study/lib/endpoint_math.R');source('study/lib/simulation_endpoints.R')
units<-read.delim(file.path(input,'units.tsv'),stringsAsFactors=FALSE)
stopifnot(nrow(units)==120L)
for(i in seq_len(nrow(units))){
 z<-readRDS(file.path(input,paste0(units$task_id[i],'.rds')))
 stopifnot(z$schema=='PUBLIC_SIMULATION_EVIDENCE_V1',z$unit$task_id==units$task_id[i])
 sim_unit_endpoints(z,out);cat('SIMULATION_REBUILD',i,'/120\n');flush.console()
}
for(p in list.files(file.path(input,'inherited'),full.names=TRUE))
 stopifnot(file.copy(p,file.path(out,basename(p)),overwrite=FALSE))
writeLines('DERIVED_ONLY_NO_MODEL_FITS',file.path(out,'DERIVED_BUILD_COMPLETE'))
