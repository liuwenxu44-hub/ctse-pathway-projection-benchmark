# Endpoint-only delta builder. Reuses frozen C1/proxy/baseline assets; no fitting.
args<-commandArgs(TRUE);stopifnot(length(args)==2);I<-normalizePath(args[1]);out<-args[2];stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
source('study/lib/cellbench_math.R')
rd<-function(n)read.delim(file.path(I,n),check.names=FALSE,stringsAsFactors=FALSE)
libs<-read.csv(file.path(I,'libraries.csv'),check.names=FALSE);el<-rd('eligibility.tsv');el$eligible<-tolower(el$eligible)=='true';tasks<-rd('canonical_tasks.tsv');plan<-rd('CONTEXT_REBUILD_PLAN.tsv');cell<-c('HCC827','H1975','H2228')
oldroot<-I
tasks$canonical_path<-file.path(I,tasks$canonical_path)
genes<-rd('reactome.tsv');reactome<-split(genes$gene_id,genes$reactome_stable_id);stopifnot(length(reactome)==1183)
c2<-rd('c2.tsv');rownames(c2)<-c2$gene_id;c2<-as.matrix(c2[,cell]);proxies<-list();truth<-list()
o<-rd('oracle.tsv');mixmap<-setNames(rep(seq_len(nrow(o)),lengths(strsplit(o$source_mixture_ids,';',fixed=TRUE))),unlist(strsplit(o$source_mixture_ids,';',fixed=TRUE)))
for(p in c('CEL','SORT')){
 c1<-read.delim(file.path(oldroot,paste0('C1_PROXY_MATERIALIZED_',p,'.tsv')),check.names=FALSE);rownames(c1)<-c1$gene_id;c1<-as.matrix(c1[,cell]);proxies[[p]]<-list(C1=c1,C2=c2)
 L<-libs[libs$platform==if(p=='CEL')'CEL-seq2' else 'SORT-seq',];W<-as.matrix(o[unname(mixmap[L$source_mix_id]),cell]);rownames(W)<-L$library_id;truth[[p]]<-W
}
scores<-list();drifts<-list();projs<-list();schema<-list();coverage<-list();scount<-0L;dcount<-0L;pcount<-0L
for(t in seq_len(nrow(tasks))){
 task<-tasks[t,];cn<-plan$context[plan$task_id==task$task_id];if(!length(cn))next
 actual<-strsplit(system2('sha256sum',shQuote(task$canonical_path),stdout=TRUE),' ')[[1]][1];stopifnot(actual==task$sha256)
 x<-readRDS(task$canonical_path);X<-x$tensor;p<-if(task$platform=='CEL-seq2')'CEL' else 'SORT';L<-libs[libs$platform==task$platform,]
 stopifnot(identical(x$libraries,L$library_id),identical(dimnames(X),list(x$genes,cell,L$library_id)),!anyDuplicated(x$genes),all(is.finite(X)))
 epic<-task$method=='EPIC-unmix';bp<-task$method=='BayesPrism'
 if(epic)stopifnot(x$scale=='conditional_log2_column_CPM13868_plus_1',x$canonical_transform=='NAMED_AXES_ONLY')
 schema[[length(schema)+1L]]<-data.frame(task=task$task_id,genes=dim(X)[1],celltypes=dim(X)[2],libraries=dim(X)[3],duplicate_genes=0,nonfinite=0,hash_match=TRUE)
 for(ctx in cn){
  U<-rd(paste0(ctx,'.tsv'))$gene_id;stopifnot(all(U%in%x$genes),!anyDuplicated(U));G<-length(U);idx<-match(U,x$genes)
  coverage[[length(coverage)+1L]]<-data.frame(task=task$task_id,context=ctx,input_genes=13868,returned_genes=length(x$genes),endpoint_genes=G,covered_reactome_genes=length(intersect(U,unique(genes$gene_id))))
  for(k in 1:3){
   M<-X[idx,k,,drop=FALSE];dim(M)<-c(G,nrow(L))
   for(pr in c('C1','C2')){
    e<-el[el$platform==task$platform & el$cell_line==cell[k] & el$proxy==pr,];stopifnot(identical(e$library,L$library_id))
    y<-proxies[[p]][[pr]][U,k];cstate<-apply(M,2,function(q)cor_reason(q,y))
    vals<-list(Spearman=as.numeric(suppressWarnings(cor(M,y,method='spearman'))),Pearson=if(epic)rep(NA_real_,nrow(L)) else as.numeric(suppressWarnings(cor(M,y,method='pearson'))))
    vals$SNE<-if(bp||epic)rep(NA_real_,nrow(L)) else sqrt(colMeans((M-y)^2))/max(sd(y),1e-8)
    for(metric in names(vals)){
     v<-vals[[metric]];st<-ifelse(is.finite(v),'AVAILABLE','UNDEFINED_CORRELATION_OR_NONFINITE')
     if(metric%in%c('Spearman','Pearson')){st<-cstate;v[cstate!='CORRELATION_DEFINED']<-NA_real_}
     if(metric=='SNE'&&bp)st[]<-'NOT_COMPARABLE_CONTRIBUTION_VS_CONDITIONAL'
     if(epic&&metric!='Spearman')st[]<-'NOT_COMPARABLE_LOG_CONDITIONAL_VS_LINEAR_CPM'
     st[!e$eligible]<-e$status[!e$eligible];v[!e$eligible]<-NA_real_
     scount<-scount+1L;scores[[scount]]<-data.frame(task=task$task_id,method=task$method,platform=task$platform,arm=task$arm,context=ctx,proxy=pr,cell_line=cell[k],composition=L$composition_block,stratum=L$stratum_id,library=L$library_id,metric,value=v,status=st,eligible=e$eligible,gene_count=G)
    }
   }
   comps<-sort(unique(L$composition_block[truth[[p]][,k]>0]));stopifnot(length(comps)==5)
   # No EPIC log-space amplitude/projection manufactured for a CPM contract.
   V<-if(!epic&&!bp)sapply(comps,function(co){ii<-which(L$composition_block==co);vecmean(M[,ii,drop=FALSE],L$stratum_id[ii])}) else NULL
   members<-lapply(reactome,function(g)match(intersect(g,U),U));counts<-lengths(members);ok<-counts>=15&counts<=500
   for(j in seq_along(comps)){
    legal<-!epic&&!bp&&!is.null(V)&&all(is.finite(V));reason<-if(epic)'NOT_COMPARABLE_LOG_CONDITIONAL_VS_LINEAR_CPM' else 'NOT_COMPARABLE_CONTRIBUTION_VS_CONDITIONAL'
    av<-NA_real_;pv<-rep(NA_real_,length(members))
    if(legal){anchor<-rowMeans(V[,-j,drop=FALSE]);delta<-V[,j]-anchor;av<-amplitude(delta,anchor);dn<-sqrt(sum(delta^2));pv[ok]<-if(dn<=1e-15)0 else vapply(members[ok],function(i)sum(delta[i])/(sqrt(length(i))*dn),0.0)}
    dcount<-dcount+1L;drifts[[dcount]]<-data.frame(task=task$task_id,method=task$method,platform=task$platform,arm=task$arm,context=ctx,cell_line=cell[k],composition=comps[j],gene_count=G,eligible_compositions=5,anchor_compositions=4,metric='DRIFT_AMPLITUDE_A',value=av,status=if(legal)'AVAILABLE' else reason)
    st<-ifelse(ok,if(legal)'AVAILABLE' else reason,'PATHWAY_SIZE_OUTSIDE_FROZEN_RANGE')
    pcount<-pcount+1L;projs[[pcount]]<-data.frame(task=task$task_id,method=task$method,platform=task$platform,arm=task$arm,context=ctx,cell_line=cell[k],composition=comps[j],pathway=names(members),member_count=counts,gene_count=G,P=pv,direction=direction(pv),threshold_exceeded=ifelse(is.finite(pv),abs(pv)>.05,NA),status=st)
   }
  }
 }
 cat('UPDATED_CONTEXTS ',task$task_id,' ',paste(cn,collapse=';'),'\n');flush.console();rm(x,X);gc(FALSE)
}
stopifnot(length(scores)>0)
write_tsv(do.call(rbind,scores),file.path(out,'CTSE_LIBRARY_LONG.tsv'))
write_tsv(do.call(rbind,drifts),file.path(out,'DRIFT_COMPOSITION_LONG.tsv'))
write_tsv(do.call(rbind,projs),file.path(out,'REACTOME_COMPOSITION_LONG.tsv'))
write_tsv(do.call(rbind,schema),file.path(out,'CANONICAL_SCHEMA_AUDIT.tsv'))
write_tsv(do.call(rbind,coverage),file.path(out,'UNIVERSE_COVERAGE.tsv'))
# Direct static proxy rows change only when their actual practical gene context changes.
dep<-rd('UNIVERSE_DEPENDENCY_AUDIT.tsv');rows<-list()
for(p in c('CEL','SORT'))for(k in 1:3){
 U<-rd(paste0(p,'_PRACTICAL.tsv'))$gene_id;x<-proxies[[p]]$C1[U,k];y<-proxies[[p]]$C2[U,k]
 rows[[length(rows)+1L]]<-data.frame(platform=if(p=='CEL')'CEL-seq2' else 'SORT-seq',context=paste0(p,'_PRACTICAL'),cell_line=cell[k],gene_count=length(U),Spearman=safe_cor(x,y),Pearson=safe_cor(x,y,'pearson'),SNE_C1_vs_C2=sne(x,y),eligible_compositions=4,interpretation='STATIC_PROXY_AGREEMENT_NOT_INDEPENDENT_REPLICATES')
}
if(length(rows))write_tsv(do.call(rbind,rows),file.path(out,'PROXY_DIRECT_AGREEMENT.tsv'))
cat('AFFECTED_ENDPOINT_DELTA_COMPLETE_NO_METHOD_FITS\n')

stopifnot(file.copy(file.path(I,'FRACTION_LIBRARY_LONG.tsv'),file.path(out,'FRACTION_LIBRARY_LONG.tsv')))
writeLines('CORRECTED_CANONICAL_DERIVED_ONLY',file.path(out,'DERIVED_BUILD_COMPLETE'))
