sim_unit_endpoints <- function(bundle, dest) {
 t <- bundle$unit; id <- as.character(t$task_id)
 ev <- bundle$evaluation; gfull <- bundle$genes; cells <- bundle$cells; samples <- bundle$samples
 grp <- as.character(ev$metadata$group)
 z <- bundle$regimes[[1]]
 old <- z$mats[c('bayesprism','unico','tca')]
 old <- old[!vapply(old,is.null,logical(1))]
 available <- lapply(old, rownames); epic <- z$mats$epic
 if (!is.null(epic)) available[['epic']] <- rownames(epic)
 stat <- data.frame(status = z$failure$epic)
 bl <- list(reference_linear=z$mats$reference, bulk_linear=z$mats$bulk,
            reference_native_log=z$log_baselines$reference, bulk_native_log=z$log_baselines$bulk)
 write_chunk <- function(rows,name) {
   if(!length(rows))return();p<-file.path(dest,name);ex<-file.exists(p)
   write.table(do.call(rbind,rows),p,sep='\t',quote=FALSE,row.names=FALSE,col.names=!ex,append=ex,na='NA')
 }
 available[['reference']]<-available[['bulk']]<-gfull
 common<-Reduce(intersect,available);common<-gfull[gfull%in%common];stopifnot(length(common)>0)
 mats<-c(old,list(epic=epic,reference=bl$reference_linear,bulk=bl$bulk_linear));methods<-c('bayesprism','unico','tca','epic','reference','bulk')
 coverage<-lapply(methods,function(m)data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,method=m,original_input_genes=2000,returned_genes=if(is.null(mats[[m]]))0 else nrow(mats[[m]]),common_genes=length(common),available=!is.null(mats[[m]]),failure=if(!is.null(mats[[m]]))'' else if(m=='epic')as.character(stat$status) else 'PRESERVED_ORIGINAL_TCA_FAILURE',celltypes=4,samples=40))
 write_chunk(coverage,'METHOD_AVAILABILITY.tsv')
 write_chunk(list(data.frame(task_id=id,gene_id=gfull,common_presence=gfull%in%common,epic_returned=gfull%in%available[['epic']])),'COMMON_RETURNED_UNIVERSE.tsv')
 rk<-list();er<-list();con<-list();dr<-list();maskrows<-list();rr<-0;ee<-0;cc<-0;dd<-0
 for(k in 1:4){
  masks<-strata_masks(ev,common,k);typechanged<-any(ev$design_changed[,k]);maskrows[[k]]<-data.frame(task_id=id,celltype=cells[k],celltype_design_changed=typechanged,all_input_changed_genes=sum(ev$design_changed[,k]),common_changed_genes=sum(masks$DESIGN_CHANGED),common_unchanged_genes=sum(masks$DESIGN_UNCHANGED))
  for(m in methods){
   x<-mats[[m]];y<-ev$truth_linear[common,k,,drop=FALSE];y<-matrix(y,nrow=length(common))
   for(s in 1:40)for(st in names(masks)){
    keep<-masks[[st]];n<-sum(keep);xx<-if(is.null(x))NULL else x[common,k,s];yy<-y[,s]
    rho<-if(is.null(xx)||n<3)NA_real_ else safe_rank(xx[keep],yy[keep])
    reason<-if(is.null(xx))'METHOD_UNAVAILABLE' else if(n<3)'INSUFFICIENT_DESIGN_STRATUM_GENES' else if(!is.finite(rho))'CONSTANT_VECTOR_RANK_UNDEFINED' else 'EVALUABLE'
    rr<-rr+1;rk[[rr]]<-data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,method=m,celltype=cells[k],celltype_design_changed=typechanged,sample_id=samples[s],group=grp[s],gene_stratum=st,genes=n,spearman=rho,status=reason)
   }
  }
  contexts<-list(LINEAR_CPM_SECONDARY=list(unico=old$unico,tca=old$tca,reference=bl$reference_linear,bulk=bl$bulk_linear),EPIC_NATIVE_LOG_SECONDARY=list(epic=epic,reference=bl$reference_native_log,bulk=bl$bulk_native_log))
  for(ctx in names(contexts)){
   truth<-if(ctx=='LINEAR_CPM_SECONDARY')ev$truth_linear else ev$truth_native_log
   Y<-matrix(truth[common,k,,drop=FALSE],nrow=length(common))
   td<-rowMeans(Y[,grp=='group1',drop=FALSE])-rowMeans(Y[,grp=='group0',drop=FALSE])
   for(m in names(contexts[[ctx]])){
    X<-contexts[[ctx]][[m]];xx<-if(is.null(X))NULL else matrix(X[common,k,,drop=FALSE],nrow=length(common))
    delta<-if(is.null(xx))NULL else rowMeans(xx[,grp=='group1',drop=FALSE])-rowMeans(xx[,grp=='group0',drop=FALSE])
    for(s in 1:40)for(st in names(masks)){
     keep<-masks[[st]];n<-sum(keep);err<-if(is.null(xx)||!n)NA_real_ else sqrt(mean((xx[keep,s]-Y[keep,s])^2));den<-if(n>=2)max(sd(Y[keep,s]),1e-8) else NA_real_
     ee<-ee+1;er[[ee]]<-data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,context=ctx,method=m,celltype=cells[k],celltype_design_changed=typechanged,sample_id=samples[s],group=grp[s],gene_stratum=st,genes=n,rmse=err,sne=err/den,status=if(is.null(xx))'METHOD_UNAVAILABLE' else if(!n)'EMPTY_DESIGN_STRATUM' else 'SAME_UNIT_LATENT_TARGET_DIAGNOSTIC')
    }
    for(st in names(masks)){
     keep<-masks[[st]];n<-sum(keep);cc<-cc+1
     con[[cc]]<-data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,context=ctx,method=m,celltype=cells[k],celltype_design_changed=typechanged,gene_stratum=st,genes=n,truth_contrast_L2=sqrt(sum(td[keep]^2)),estimated_contrast_L2=if(is.null(delta)||!n)NA_real_ else sqrt(sum(delta[keep]^2)),contrast_error_L2=if(is.null(delta)||!n)NA_real_ else sqrt(sum((delta[keep]-td[keep])^2)))
    }
    for(path in c('signal','negative_control')){
     p<-if(is.null(delta))c(P=NA_real_,norm=NA_real_,members=sum(common%in%sprintf('sg%04d',if(path=='signal')1:50 else 51:100))) else projection(delta,common,path);tp<-projection(td,common,path)
     dd<-dd+1;dr[[dd]]<-data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,context=ctx,method=m,celltype=cells[k],celltype_design_changed=typechanged,pathway=path,genes=length(common),pathway_members=unname(p['members']),estimated_P=unname(p['P']),truth_P=unname(tp['P']),estimated_direction=unname(classify(p['P'])),truth_direction=unname(classify(tp['P'])),evaluable=is.finite(p['P']),semantic_role='SCALE_SPECIFIC_SECONDARY_NOT_FOUR_METHOD_RANKING')
    }
   }
   for(path in c('signal','negative_control')){
    tp<-projection(td,common,path);dd<-dd+1;dr[[dd]]<-data.frame(task_id=id,scenario=t$scenario,replicate_id=t$replicate_id,context=ctx,method='always_zero',celltype=cells[k],celltype_design_changed=typechanged,pathway=path,genes=length(common),pathway_members=unname(tp['members']),estimated_P=0,truth_P=unname(tp['P']),estimated_direction='zero',truth_direction=unname(classify(tp['P'])),evaluable=TRUE,semantic_role='ANALYTIC_ZERO_PREDICTION_NO_MODEL')
   }
  }
 }
 write_chunk(rk,'GENE_RANK_RECOVERY_LONG.tsv');write_chunk(er,'SAME_SCALE_ERROR_LONG.tsv');write_chunk(con,'SCALE_SPECIFIC_CONTRAST_LONG.tsv');write_chunk(dr,'SCALE_SPECIFIC_DIRECTION_LONG.tsv');write_chunk(maskrows,'DESIGN_CHANGE_STRATA.tsv')
}
