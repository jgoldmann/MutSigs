
fit_to_signatures <- function(mut_matrix, signatures) {
  require(nnls)
  if(!(dim(mut_matrix)[1] == 96))stop("Mutation count matrix input should have dimensions 96 X n samples")
  if(any(is.null(colnames(mut_matrix))))stop("Mutation count matrix must have proper rownames")
  if(!(dim(signatures)[1] == 96))stop("Signatures input should have dimensions 96 X n signatures")
  lsq_fits <- apply(mut_matrix, 2, function(x){nnls(signatures, x)})

  res <-
    tibble(
      contrib = foreach(j=1:length(lsq_fits), .combine = c) %do% {lsq_fits[[j]]$x},
      group = foreach(samplename = colnames(mut_matrix), .combine = c) %do% {rep(samplename, ncol(signatures))},
      groupMutations = foreach(j = 1:ncol(mut_matrix), .combine = c) %do% {rep(sum(mut_matrix[,j]), ncol(signatures))},
      signature = foreach(samplename = colnames(mut_matrix), .combine = c) %do% {colnames(signatures)},
      coefDeter = foreach(j=1:length(lsq_fits), .combine = c) %do% {1-(sum(lsq_fits[[j]]$residuals^2)/sum(lsq_fits[[j]]$fitted^2)) %>% rep(ncol(signatures))},
      residualGini = foreach(j=1:length(lsq_fits), .combine = c) %do% {
        lsq_fits[[j]]$residuals %>%
          abs() %>%
          edgeR::gini() %>%
          rep(ncol(signatures))}
    )

  res <-
    bind_rows(res,
              tibble(
                contrib = foreach(j=1:length(lsq_fits), .combine = c) %do% {sum(lsq_fits[[j]]$residuals)},
                group = colnames(mut_matrix),
                groupMutations = foreach(j = 1:ncol(mut_matrix), .combine = c) %do% {sum(mut_matrix[,j])},
                signature = "unfittable",
                coefDeter = foreach(j=1:length(lsq_fits), .combine = c) %do% {1-(sum(lsq_fits[[j]]$residuals^2)/sum(lsq_fits[[j]]$fitted^2))},
                residualGini = foreach(j=1:length(lsq_fits), .combine = c) %do% {
                  lsq_fits[[j]]$residuals %>%
                    abs() %>%
                    edgeR::gini()}
              )) %>%
    arrange(group)

  return(res)
}


resampleContributions <- function(phasedDNMs, signatures, n, group, rnd.seed = 1234) {
  stopifnot(is.data.frame(phasedDNMs) & all(c("surrounding", "substitution") %in% colnames(phasedDNMs)))

  contributions <-
    foreach(i = 1:n,
            .combine = bind_rows) %do% {
      set.seed(rnd.seed*i)
      DNMsubs <- phasedDNMs[sample(1:nrow(phasedDNMs), replace = TRUE), , drop=FALSE]

      triNucSubs <- as.data.frame(table(DNMsubs$surrounding,
                                        DNMsubs$substitution,
                                        DNMsubs %>% pull(!!group)))
      colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

      triNucSubs %>%
        unite(mutation, substitution, context) %>%
        spread_(group, "frequency") %>%
        dplyr::select(-mutation) %>%
        as.matrix() %>%
        fit_to_signatures(signatures) %>%
        mutate(i = i)
            }

  contributions %>%
    group_by(group, signature) %>%
    summarise(ci95upper = quantile(contrib, 0.975),
              ci95lower = quantile(contrib, 0.025),
              sdv = sd(contrib)) %>%
    return()
}


#' Main function for fitting signatures by group with confidence intervals.
#'
#' @param DNMs A data.frame that contains information about the mutations.
#'        Necessary columns are \code{surrounding}, \code{substitution} plus a column  with information about the groups.
#' @param signatures A matrix with relative abundances of all 96 mutation types in the rows and all signatures to be investigated in the columns.
#' @param group Name of the column in \code{DNMs} to base the groups upon.
#' @param rnd.seed Random seed to involve in the resampling of the fitting.
#' @param n Number of resamplings to be done.
assessSignatures <- function(DNMs,
                              signatures = cancer_signatures,
                              group = "parent",
                              rnd.seed = 1234,
                              n=1000) {
  require(tidyverse)
  require(nnls)
  require(foreach)

  stopifnot(c("surrounding", "substitution", group) %in% colnames(DNMs))
  stopifnot(dim(signatures)[[1]]==96)
  stopifnot(is.numeric(rnd.seed))
  stopifnot(is.numeric(n))

  triNucSubs <- as.data.frame(table(DNMs$surrounding,
                                    DNMs$substitution,
                                    DNMs %>% pull(!!group)))
  colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

  resampleContributions(DNMs, signatures, 100, group) %>%
    full_join(
      triNucSubs %>%
        unite(mutation, substitution, context) %>%
        spread_(group, "frequency") %>%
        dplyr::select(-mutation) %>%
        as.matrix() %>%
        fit_to_signatures(signatures),
      by = c("group", "signature")) %>%
    mutate(relC = contrib/groupMutations,
           relci95upper = ci95upper/groupMutations,
           relci95lower = ci95lower/groupMutations,
           relSd = sdv/groupMutations) %>%
    return()
}





