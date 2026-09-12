# Typed scale helpers: no inference from numeric magnitude and no model calls.
ctse_epic_linear_input <- function(x, gene_ids, sample_ids) {
  ctse_assert(!inherits(x, "ctse_epic_log_response"), "DOUBLE_LOG_TRANSFORM_FORBIDDEN")
  ctse_assert(is.matrix(x) && is.numeric(x) && all(is.finite(x)) && all(x>=0), "LINEAR_INPUT_INVALID")
  ctse_assert(identical(dimnames(x),list(gene_ids,sample_ids)) && !anyDuplicated(gene_ids) &&
                !anyDuplicated(sample_ids) && all(colSums(x)>0), "FULL_INPUT_AXIS_INVALID")
  structure(list(values=x, full_gene_ids=gene_ids, sample_ids=sample_ids,
                 scale="DECLARED_LINEAR_FULL_INPUT"),class="ctse_epic_linear_input")
}
ctse_epic_log_response <- function(x) {
  ctse_assert(inherits(x,"ctse_epic_linear_input"), "DECLARED_LINEAR_INPUT_REQUIRED_NO_DOUBLE_LOG")
  den<-colSums(x$values)
  L<-log2(sweep(x$values,2,den,"/")*1e6+1)
  structure(list(values=L,full_gene_ids=x$full_gene_ids,sample_ids=x$sample_ids,
                 full_denominator=den,transform_count=1L,scale="log2_full_support_CPM_plus_1"),
            class="ctse_epic_log_response")
}
ctse_epic_subset_response <- function(x, returned_gene_ids) {
  ctse_assert(inherits(x,"ctse_epic_log_response") && x$transform_count==1L,"LOG_CONTRACT_REQUIRED")
  ctse_assert(!anyDuplicated(returned_gene_ids) && length(returned_gene_ids)>0 &&
                all(returned_gene_ids%in%rownames(x$values)),"RETURNED_UNIVERSE_INVALID")
  x$values<-x$values[returned_gene_ids,,drop=FALSE]
  x # Full denominator and full-gene identity intentionally unchanged.
}
ctse_validate_fit_access <- function(input, arm) {
  ctse_assert(arm%in%c("practical","oracle"),"FRACTION_ARM_REQUIRED")
  scan<-function(z){
    if(!is.list(z))return(character())
    c(names(z),unlist(lapply(z,scan),use.names=FALSE))
  }
  fields<-tolower(scan(input))
  forbidden<-c("truth","evaluation","c1_proxy","c2_proxy","performance","endpoint")
  if(arm=="practical")forbidden<-c(forbidden,"oracle","known_proportions")
  ctse_assert(!any(vapply(forbidden,function(p)any(grepl(p,fields,fixed=TRUE)),logical(1))),
              "EVALUATION_OR_ORACLE_ACCESS_FORBIDDEN")
  invisible(TRUE)
}
