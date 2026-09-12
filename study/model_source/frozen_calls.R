# Source-only public interfaces to the recorded calls. No call runs on source().
# They are NOT invoked by repro/run.py or installed in the derived-only image.
# See CALL_SOURCE_MAP.tsv and docs/MODEL_EXECUTION.md before any new fit.
study_fit_gate <- function(method, input, arm, allow_new_fit) {
  stopifnot(isTRUE(allow_new_fit))
  ctse_validate_fit_access(input, arm)
  expected <- c(BayesPrism="2.2.3",Unico="0.1.0",TCA="1.2.1",EPICunmix="0.0.1")
  stopifnot(method %in% names(expected),
            as.character(utils::packageVersion(method))==expected[[method]])
  stopifnot(identical(rownames(input$CPM),input$genes),
            identical(colnames(input$CPM),input$sample_ids),
            identical(dimnames(input$W),list(input$sample_ids,input$cell_ids)))
}
study_rng <- function(seed) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(as.integer(seed))
}
simulation_unico_call <- function(input, seed, allow_new_fit=FALSE) {
  study_fit_gate("Unico",input,"practical",allow_new_fit);study_rng(seed)
  fit <- Unico::Unico(input$CPM,input$W,C1=NULL,C2=NULL,parallel=FALSE,
                      log_file=NULL,verbose=FALSE)
  native <- Unico::tensor(input$CPM,input$W,C1=NULL,C2=NULL,Unico.mdl=fit,
                         parallel=FALSE,log_file=NULL,verbose=FALSE)
  list(raw_fit=fit,raw_output=native,
       tensor=q2_canonicalize_unico(native,input$genes,input$cell_ids,input$sample_ids),
       scale="signed_conditional_CPM")
}
simulation_tca_call <- function(input, seed, allow_new_fit=FALSE) {
  study_fit_gate("TCA",input,"practical",allow_new_fit);study_rng(seed)
  fit <- TCA::tca(input$CPM,input$W,C1=NULL,C2=NULL,refit_W=FALSE,parallel=FALSE,
                  max_iters=10L,log_file=NULL,verbose=FALSE)
  native <- TCA::tensor(input$CPM,tca.mdl=fit,scale=FALSE,parallel=FALSE,
                       log_file=NULL,verbose=FALSE)
  list(raw_fit=fit,raw_output=native,
       tensor=q2_canonicalize_tca(native,input$genes,input$cell_ids,input$sample_ids),
       scale="signed_conditional_CPM")
}
simulation_bayesprism_call <- function(input,reference,seed,allow_new_fit=FALSE) {
  study_fit_gate("BayesPrism",input,"practical",allow_new_fit)
  ctse_validate_fit_access(reference,"practical");study_rng(seed)
  stopifnot(identical(rownames(reference$counts),input$genes),
            identical(colnames(reference$counts),as.character(reference$metadata$cell_id)))
  object <- BayesPrism::new.prism(reference=t(reference$counts),input.type="count.matrix",
    cell.type.labels=as.character(reference$metadata$cell_type),
    cell.state.labels=paste(reference$metadata$cell_type,reference$metadata$state,sep="__"),
    key=NULL,mixture=t(input$counts))
  fit <- BayesPrism::run.prism(object,n.cores=1L)
  expression <- array(NA_real_,c(length(input$genes),length(input$cell_ids),length(input$sample_ids)),
    dimnames=list(input$genes,input$cell_ids,input$sample_ids))
  for(cell in input$cell_ids)expression[,cell,] <- t(bp_coerce_public_expression_matrix(
    BayesPrism::get.exp(fit,"type",cell),input$sample_ids,input$genes))
  normalized <- bp_normalize_type_expression(expression,input$sample_ids,input$cell_ids,input$genes)
  list(raw_fit=fit,expression_raw=expression,tensor=normalized$normalized,
       raw_sums=normalized$raw_sums,scale="type_profile_unit_sum_probability")
}
cellbench_unico_call <- function(input,arm,seed,log_directory,allow_new_fit=FALSE) {
  study_fit_gate("Unico",input,arm,allow_new_fit);study_rng(seed)
  stopifnot(length(input$genes)==13868L,"ENSG00000105519"%in%input$genes,
            dir.exists(log_directory))
  execution_genes <- setdiff(input$genes,"ENSG00000105519")
  X <- input$CPM[execution_genes,,drop=FALSE];W <- input$W
  fit <- Unico::Unico(X=X,W=W,C1=NULL,C2=NULL,fit_tau=FALSE,mean_penalty=0,
    var_penalty=0.1,covar_penalty=0.1,mean_max_iterations=2,var_max_iterations=3,
    nloptr_opts_algorithm="NLOPT_LN_COBYLA",init_weight="default",max_u=1,max_v=1,
    parallel=FALSE,num_cores=1,log_file=file.path(log_directory,"Unico.log"),verbose=TRUE,debug=FALSE)
  native <- Unico::tensor(X=X,W=W,C1=NULL,C2=NULL,Unico.mdl=fit,parallel=FALSE,num_cores=1,
    log_file=file.path(log_directory,"Unico_tensor.log"),verbose=TRUE,debug=FALSE)
  tensor <- aperm(native,c(2L,1L,3L))[execution_genes,input$cell_ids,input$sample_ids,drop=FALSE]
  list(raw_fit=fit,raw_output=native,tensor=tensor,scale="signed_native_CPM")
}
cellbench_tca_call <- function(input,arm,seed,log_directory,allow_new_fit=FALSE) {
  study_fit_gate("TCA",input,arm,allow_new_fit);study_rng(seed);stopifnot(dir.exists(log_directory))
  fit <- TCA::tca(X=input$CPM,W=input$W,C1=NULL,C2=NULL,refit_W=FALSE,tau=NULL,
    vars.mle=FALSE,constrain_mu=FALSE,parallel=FALSE,num_cores=1,max_iters=10,
    log_file=file.path(log_directory,"TCA.log"),verbose=TRUE)
  native <- TCA::tensor(X=input$CPM,tca.mdl=fit,scale=FALSE,parallel=FALSE,num_cores=1,
    log_file=file.path(log_directory,"TCA_tensor.log"),verbose=TRUE)
  stopifnot(length(native)==3L)
  names(native) <- if(is.null(names(native)))input$cell_ids else names(native)
  tensor <- array(NA_real_,c(length(input$genes),3L,length(input$sample_ids)),
    dimnames=list(input$genes,input$cell_ids,input$sample_ids))
  for(cell in input$cell_ids)tensor[,cell,] <- native[[cell]][input$genes,input$sample_ids,drop=FALSE]
  list(raw_fit=fit,raw_output=native,tensor=tensor,scale="signed_native_CPM")
}
cellbench_epic_corrected_stage2_call <- function(a,initial,seed,arm,allow_new_fit=FALSE) {
  stopifnot(isTRUE(allow_new_fit),as.character(utils::packageVersion("EPICunmix"))=="0.0.1")
  ctse_validate_fit_access(a,arm);ctse_validate_fit_access(initial,arm)
  stopifnot(identical(dimnames(a$bulk),list(a$genes,a$libraries)),
    identical(dimnames(initial$A),list(a$genes,a$cell_lines,a$libraries)),
    identical(dimnames(a$W),list(a$libraries,a$cell_lines)),length(a$genes)==13868L,
    all(is.finite(a$bulk)),all(is.finite(a$W)))
  # a$bulk is the separately persisted, scale-verified response. No transformation here.
  install_shim();set.seed(as.integer(seed))
  result <- EPICunmix::run_epic_unmix(bulk=a$bulk,frac=a$W,input_cts=initial,
    outf=FALSE,nstop=1,delta=.1,nu0=50,nu1=50,seed=as.integer(seed),ncore=1)
  X<-result$A;genes<-a$genes[a$genes%in%rownames(X)]
  X<-X[genes,a$cell_lines,a$libraries,drop=FALSE]
  stopifnot(all(is.finite(X)),!anyDuplicated(genes))
  list(raw_output=result,tensor=X,scale="conditional_log2_column_CPM13868_plus_1")
}
