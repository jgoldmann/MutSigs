

resampleContributions <- function(phasedDNMs, cancer_signatures, n, group = "parent", rnd.seed = 1234) {
  stopifnot(is.numeric(n) & length(n) == 1)
  stopifnot(is.data.frame(phasedDNMs) & all(c("surrounding", "substitution", group) %in% colnames(phasedDNMs)))

  sampleAndFit <- function(i){
    set.seed(rnd.seed*i)
    DNMsubs <- phasedDNMs[sample(1:nrow(phasedDNMs), replace = TRUE), , drop=FALSE]
    DNMsubs_tbl <- as.data.frame(table("context" = DNMsubs$surrounding,
                                       "substitution" = DNMsubs$substitution,
                                       "group" = DNMsubs[,group]))
    fit_res_i <-
      DNMsubs_tbl %>%
      unite(mutation, substitution, context) %>%
      spread(group, Freq) %>%
      dplyr::select(-mutation) %>%
      fit_to_signatures(cancer_signatures)

    return(foreach(j=1:length(fit_res_i), .combine = cbind) %do% {fit_res_i[[j]]$x})
  }

  contributions <-
    foreach(i=1:n) %do% sampleAndFit(i)

  ctbA <- abind(contributions, along=3)
  return(ctbA)
}




fit_to_signatures <- function(mut_matrix, signatures) {
  require(nnls)
  if(!(dim(mut_matrix)[1] == 96))stop("Mutation count matrix input should have dimensions 96 X n samples")
  if(!(dim(mut_matrix)[1] == 96))stop("Signatures input should have dimensions 96 X n signatures")
  lsq_fits <- apply(mut_matrix, 2, function(x){nnls(signatures, x)})
  return(lsq_fits)
}




#' Main function for fitting signatures by group with confidence intervals.
#'
#' @param DNMs A data.frame that contains information about the mutations.
#'        Necessary columns are \code{surrounding}, \code{substitution} plus a column  with information about the groups.
#' @param cancer_signatures A matrix with relative abundances of all 96 mutation types in the rows and all signatures to be investigated in the columns.
#' @param group Name of the column in \code{DNMs} to base the groups upon.
#' @param rnd.seed Random seed to involve in the resampling of the fitting.
#' @param n Number of resamplings to be done.
assessSignatures <- function(DNMs,
                             cancer_signatures = cancer_signatures,
                             group = "parent",
                             rnd.seed = 1234,
                             n=1000) {
  require(tidyverse)
  require(nnls)
  require(foreach)
  stopifnot(c("surrounding", "substitution", group) %in% colnames(DNMs))
  stopifnot(dim(cancer_signatures)[[1]]==96)
  stopifnot(is.integer(rnd.seed))
  stopifnot(is.infinite(n))

  triNucSubs <- as.data.frame(table(DNMs$surrounding,
                                    DNMs$substitution,
                                    DNMs[,group]))
  colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

  fit_res <-
    triNucSubs %>%
    unite(mutation, substitution, context) %>%
    spread_(group, "frequency") %>%
    dplyr::select(-mutation) %>%
    fit_to_signatures(cancer_signatures)
  fit_resContribs <- foreach(j=1:length(fit_res), .combine = cbind) %do% {fit_res[[j]]$x}
  colnames(fit_resContribs) <- names(fit_res)

  contributionsArray <- resampleContributions(DNMs, cancer_signatures, n, group, rnd.seed = rnd.seed)
  contributionsSd <- apply(contributionsArray, c(1,2), sd)
  contributions95lb <- apply(contributionsArray, c(1,2), function(x){quantile(x, 0.975)})
  contributions95ub <- apply(contributionsArray, c(1,2), function(x){quantile(x, 0.025)})

  ctrDf <-
    cbind(gather_(as.data.frame(fit_resContribs),  group, "contrib",  colnames(fit_resContribs)),
          gather_(as.data.frame(contributionsSd), group, "sdv",      colnames(contributionsSd)),
          gather_(as.data.frame(contributions95lb), group, "ci95upper",      colnames(contributions95lb)),
          gather_(as.data.frame(contributions95ub), group, "ci95lower",      colnames(contributions95ub)))
  ctrDf$signature <- colnames(cancer_signatures)

  relContribs <-
    ctrDf %>%
    setNames(.,make.unique(names(.))) %>%
    group_by_(group) %>%
    dplyr::mutate(relC=contrib/sum(contrib),
                  relSd=sdv/sum(contrib))
  relContribs$signature <-  factor(relContribs$signature,
                                   levels = unique(relContribs$signature)[order(as.numeric(substring(unique(relContribs$signature), 11)))],
                                   ordered = TRUE)
  relDeviance <- sapply(fit_res, function(x){x$deviance})/table(DNMs[,group])
  relContribs <- full_join(data.frame(relDeviance=as.vector(relDeviance),
                                      group=names(relDeviance)),
                           relContribs,
                           by=c(group=group))
  relContribs <-
    relContribs %>%
    group_by(group) %>%
    mutate(relci95upper = ci95upper/sum(contrib),
           relci95lower = ci95lower/sum(contrib)) %>%
    ungroup()
  colnames(relContribs)[colnames(relContribs)=="group"] <- group
  return(relContribs)
}




