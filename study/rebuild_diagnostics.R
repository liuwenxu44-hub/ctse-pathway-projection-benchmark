# Rebuild from saved draws. CTS is reconstructed exactly from Sol, not resampled.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
input<-normalizePath(args[1]);out<-args[2];stopifnot(!file.exists(out));dir.create(out,recursive=TRUE)
dir.create(file.path(out,'JOB_DIAGNOSTICS'))
source('study/lib/epic_diagnostics.R');source('study/lib/cts_from_draws.R')
stopifnot(as.character(packageVersion('posterior'))=='1.7.0')
rd<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
sha<-function(p)unname(tools::sha256sum(p))
tsv<-function(x,p)write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
append_rows<-function(x,name){p<-file.path(out,name);ex<-file.exists(p);write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,col.names=!ex,append=ex,na='NA')}
jobs<-rd(file.path(input,'jobs.tsv'));genes<-rd(file.path(input,'genes.tsv'))$gene_id
stopifnot(nrow(jobs)==24L,length(genes)==20L)
inventory<-list();checks<-list()
for(ji in seq_len(nrow(jobs))){
 j<-jobs[ji,];jd<-file.path(input,j$job_id);rows<-list()
 for(gi in seq_along(genes)){
  gene<-genes[gi];records<-list();descriptive<-list()
  for(chain in 1:4){
   p<-file.path(jd,paste0(gene,'__chain',chain,'.rds'));z<-readRDS(p)
   expected_seed<-as.integer(300000000+100000*j$unit_index+10000*j$stage+1000*j$schedule_index+10*gi+chain)
   stopifnot(z$effective_seed==expected_seed,z$gene==gene,z$stage==j$stage,z$actual_sampler_starts==1L)
   z$CTS<-epic_cts_draws(z$Sol,z$sample_ids,z$cell_ids);records[[chain]]<-z
   inventory[[length(inventory)+1L]]<-data.frame(job_id=j$job_id,chain_id=paste(j$job_id,gene,paste0('chain',chain),sep='__'),gene_id=gene,chain=chain,status=z$original_status,fit_count=z$actual_sampler_starts,raw_sha256=z$original_raw_sha256,effective_seed=z$effective_seed)
  }
  z<-epic_diagnose_chains(records);z$gene_id<-gene;z$job_id<-j$job_id
  for(q in c('Sol','VCV','CTS')){
   a<-do.call(rbind,lapply(records,function(r)r[[q]]));zz<-z[z$quantity==q,]
   stopifnot(identical(zz$parameter,colnames(a)))
   descriptive[[length(descriptive)+1L]]<-data.frame(unit_id=j$unit_id,unit_index=j$unit_index,stage=j$stage,schedule=j$schedule,job_id=j$job_id,gene_id=gene,quantity=q,parameter=colnames(a),posterior_mean=colMeans(a),posterior_sd=zz$posterior_sd,mcse_mean=zz$mcse_mean,rhat=zz$rhat,ess_bulk=zz$ess_bulk,ess_tail=zz$ess_tail,status=zz$status,row.names=NULL)
  }
  append_rows(do.call(rbind,descriptive),'POSTERIOR_SUMMARIES_LONG.tsv')
  rows[[gi]]<-z;rm(records);gc(FALSE)
 }
 regenerated<-file.path(out,'JOB_DIAGNOSTICS',paste0(j$job_id,'.tsv'));tsv(do.call(rbind,rows),regenerated)
 persisted<-file.path(jd,'expected_diagnostics.tsv');exact<-sha(regenerated)==sha(persisted)
 checks[[ji]]<-data.frame(job_id=j$job_id,regenerated_sha256=sha(regenerated),original_sha256=sha(persisted),status=if(exact)'BYTE_IDENTICAL'else'PLATFORM_DIFFERENCE_REQUIRES_VALUE_AUDIT')
 a<-do.call(rbind,rows);a$unit_id<-j$unit_id;a$stage<-j$stage;a$schedule<-j$schedule
 append_rows(a,'ALL_PARAMETER_DIAGNOSTICS.tsv');cat('DIAGNOSTIC_REBUILD',ji,'/24\n');flush.console()
}
tsv(do.call(rbind,inventory),file.path(out,'CHAIN_EXECUTION_INVENTORY.tsv'))
tsv(do.call(rbind,checks),file.path(out,'PERSISTED_DIAGNOSTIC_REBUILD_CHECK.tsv'))
writeLines('DERIVED_FROM_FROZEN_CHAINS_NO_SAMPLER_CALLS',file.path(out,'DERIVED_BUILD_COMPLETE'))
