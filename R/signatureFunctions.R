
fit_to_signatures2 <- function(mut_matrix, signatures) {
  require(nnls)
  if(!(dim(mut_matrix)[1] == 96))stop("Mutation count matrix input should have dimensions 96 X n samples")
  if(!(dim(mut_matrix)[1] == 96))stop("Signatures input should have dimensions 96 X n signatures")
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


resampleContributions2 <- function(phasedDNMs, signatures, n, group, rnd.seed = 1234) {
  stopifnot(is.data.frame(phasedDNMs) & all(c("surrounding", "substitution") %in% colnames(phasedDNMs)))

  contributions <-
    foreach(i = 1:n,
            .combine = bind_rows) %do% {
      set.seed(rnd.seed*i)
      DNMsubs <- phasedDNMs[sample(1:nrow(phasedDNMs), replace = TRUE), , drop=FALSE]

      triNucSubs <- as.data.frame(table(DNMsubs$surrounding,
                                        DNMsubs$substitution,
                                        DNMsubs %>% pull(group)))
      colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

      triNucSubs %>%
        unite(mutation, substitution, context) %>%
        spread_(group, "frequency") %>%
        dplyr::select(-mutation) %>%
        fit_to_signatures2(signatures) %>%
        mutate(i = i)
            }

  contributions %>%
    group_by(group, signature) %>%
    summarise(ci95upper = quantile(contrib, 0.975),
              ci95lower = quantile(contrib, 0.025),
              sdv = sd(contrib)) %>%
    return()
}

assessSignatures2 <- function(DNMs,
                              signatures = cancer_signatures,
                              group = "parent",
                              rnd.seed = 1234,
                              n=1000) {
  require(tidyverse)
  require(nnls)
  require(foreach)
  require(abind)

  stopifnot(c("surrounding", "substitution", group) %in% colnames(DNMs))
  stopifnot(dim(signatures)[[1]]==96)
  stopifnot(is.numeric(rnd.seed))
  stopifnot(is.numeric(n))

  triNucSubs <- as.data.frame(table(DNMs$surrounding,
                                    DNMs$substitution,
                                    DNMs %>% pull(group)))
  colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

  resampleContributions2(phasedDNMs, signatures, 100, group = "group") %>%
    full_join(
      triNucSubs %>%
        unite(mutation, substitution, context) %>%
        spread_(group, "frequency") %>%
        dplyr::select(-mutation) %>%
        fit_to_signatures2(signatures),
      by = c("group", "signature")) %>%
    mutate(relC = contrib/groupMutations,
           relci95upper = ci95upper/groupMutations,
           relci95lower = ci95lower/groupMutations,
           relSd = sdv/groupMutations) %>%
    return()
}





resampleContributions <- function(phasedDNMs, cancer_signatures, n, group = "parent", rnd.seed = 1234) {
  stopifnot(is.numeric(n) & length(n) == 1)
  if(!is.na(group)){
    stopifnot(is.data.frame(phasedDNMs) & all(c("surrounding", "substitution", group) %in% colnames(phasedDNMs)))
  } else {
    stopifnot(is.data.frame(phasedDNMs) & all(c("surrounding", "substitution") %in% colnames(phasedDNMs)))
  }

  sampleAndFit <- function(i){
    set.seed(rnd.seed*i)
    DNMsubs <- phasedDNMs[sample(1:nrow(phasedDNMs), replace = TRUE), , drop=FALSE]
    if(is.na(group)) {
      DNMsubs_tbl <- as.data.frame(table("context" = DNMsubs$surrounding,
                                         "substitution" = DNMsubs$substitution))
      fit_res_i <-
        DNMsubs_tbl %>%
        unite(mutation, substitution, context) %>%
        dplyr::select(-mutation) %>%
        fit_to_signatures(cancer_signatures)
    } else {
      DNMsubs_tbl <- as.data.frame(table("context" = DNMsubs$surrounding,
                                         "substitution" = DNMsubs$substitution,
                                         "group" = DNMsubs %>% pull(group)))
      fit_res_i <-
        DNMsubs_tbl %>%
        unite(mutation, substitution, context) %>%
        spread(group, Freq) %>%
        dplyr::select(-mutation) %>%
        fit_to_signatures(cancer_signatures)
    }
    return(foreach(j=1:length(fit_res_i), .combine = cbind) %do% {fit_res_i[[j]]$x})
  }

  contributions <-
    foreach(i=1:n) %do% sampleAndFit(i)

  dims <- if_else(is.na(group), 2, 3)
  ctbA <- abind(contributions, along=dims)
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
  require(abind)

  stopifnot(c("surrounding", "substitution", group) %in% colnames(DNMs))
  stopifnot(dim(signatures)[[1]]==96)
  stopifnot(is.numeric(rnd.seed))
  stopifnot(is.numeric(n))

  triNucSubs <- as.data.frame(table(DNMs$surrounding,
                                    DNMs$substitution,
                                    DNMs %>% pull(group)))
  colnames(triNucSubs) <- c("context", "substitution", group, "frequency")

  fit_res <-
    triNucSubs %>%
    unite(mutation, substitution, context) %>%
    spread_(group, "frequency") %>%
    dplyr::select(-mutation) %>%
    fit_to_signatures(signatures)
  fit_resContribs <- foreach(j=1:length(fit_res), .combine = cbind) %do% {fit_res[[j]]$x}
  colnames(fit_resContribs) <- names(fit_res)

  contributionsArray <- resampleContributions(DNMs, signatures, n, group, rnd.seed = rnd.seed)
  contributionsSd <- apply(contributionsArray, c(1,2), sd)
  contributions95lb <- apply(contributionsArray, c(1,2), function(x){quantile(x, 0.975)})
  contributions95ub <- apply(contributionsArray, c(1,2), function(x){quantile(x, 0.025)})

  ctrDf <-
    cbind(gather_(as.data.frame(fit_resContribs),  group, "contrib",  colnames(fit_resContribs)),
          gather_(as.data.frame(contributionsSd), group, "sdv",      colnames(contributionsSd)),
          gather_(as.data.frame(contributions95lb), group, "ci95upper",      colnames(contributions95lb)),
          gather_(as.data.frame(contributions95ub), group, "ci95lower",      colnames(contributions95ub)))
  ctrDf$signature <- colnames(signatures)

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







assessSingleGroup <- function(testData,
                  signatures = cancer_signatures,
                  rnd.seed = 1234,
                  n=100) {
  require(tidyverse)
  require(nnls)
  require(foreach)
  require(abind)

  stopifnot(c("surrounding", "substitution") %in% colnames(testData))
  stopifnot(dim(signatures)[[1]]==96)
  stopifnot(is.numeric(rnd.seed))
  stopifnot(is.numeric(n))

  triNucSubs <- as.data.frame(table(testData$surrounding,
                                    testData$substitution))
  colnames(triNucSubs) <- c("context", "substitution", "frequency")
  fit_res <-
    triNucSubs %>%
    unite(mutation, substitution, context) %>%
    dplyr::select(-mutation) %>%
    fit_to_signatures(signatures)
  fit_resContribs <- foreach(j=1:length(fit_res), .combine = cbind) %do% {fit_res[[j]]$x}

  contributionsArray <- resampleContributions(testData, signatures, n, group=NA, rnd.seed = rnd.seed)
  contributionsSd <- apply(contributionsArray, 1, sd)
  contributions95lb <- apply(contributionsArray, 1, function(x){quantile(x, 0.975)})
  contributions95ub <- apply(contributionsArray, 1, function(x){quantile(x, 0.025)})

  ctrDf <-
    data.frame("contrib" = fit_resContribs,
               "sdv" = contributionsSd,
               "ci95upper" = contributions95lb,
               "ci95lower" = contributions95ub)
  ctrDf$signature <- colnames(signatures)
  relContribs <-
    ctrDf %>%
    dplyr::mutate(relC=contrib/sum(contrib),
                  relSd=sdv/sum(contrib))

  relContribs$signature <-  factor(relContribs$signature,
                                   levels = unique(relContribs$signature)[order(as.numeric(substring(unique(relContribs$signature), 11)))],
                                   ordered = TRUE)
  relContribs %>%
    mutate(relci95upper = ci95upper/sum(contrib),
           relci95lower = ci95lower/sum(contrib),
           relDeviance = fit_res$frequency$deviance/nrow(testData)) %>%
    return()
}


