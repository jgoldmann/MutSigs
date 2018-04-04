library(tidyverse)
library(MutSigs)

context("assesSignatures")

possibleSurroundings <-
  expand.grid(
    fistLetter=c("A", "C", "G", "T"),
    secondLetter=".",
    thirdLetter=c("A", "C", "G", "T")) %>%
  unite("surrounding", fistLetter, secondLetter, thirdLetter, sep="") %>%
  pull(surrounding)

possibleSubstitutions <-
  c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

n <- 1000
set.seed(1234)

testData <- data.frame(
  surrounding = sample(possibleSurroundings, n, replace = TRUE),
  substitution = sample(possibleSubstitutions, n, replace = TRUE),
  group = sample(c("A", "B"),                n, replace = TRUE)
)

result <- suppressWarnings(assessSignatures(testData,
                                            group="group"))

test_that("the results are ok", {
  expect_true(all(dim(result)==c(60, 14)))
  expect_true(all(c("relDeviance", "group", "contrib", "sdv", "relC", "relSd", "relci95upper", "relci95lower") %in% colnames(result)))
  expect_equal(sum(result$relC),2)
  expect_equal(head(result$contrib,4), c(37.3, 0., 242., 0.), tolerance = 1)
  expect_equal(head(result$sdv,4), c(13.7, 0.561, 30.5, 7.24), tolerance = 0.5)
})
