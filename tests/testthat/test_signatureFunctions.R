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

test_that("grouped sets are derived ok", {
  expect_true(all(dim(result)==c(60, 14)))
  expect_true(all(c("relDeviance", "group", "contrib", "sdv", "relC", "relSd", "relci95upper", "relci95lower") %in% colnames(result)))
  expect_equal(sum(result$relC),2)
  expect_equal(head(result$contrib,4), c(37.3, 0., 242., 0.), tolerance = 1)
  expect_equal(head(result$sdv,4), c(13.7, 0.561, 30.5, 7.24), tolerance = 0.5)
})

resultSingle <- assessSingleGroup(testData %>% filter(group=="A"))

test_that("single sets are derived ok", {
  expect_true(all(dim(resultSingle)==c(30, 10)))
  expect_true(all(c("relDeviance", "contrib", "sdv", "relC", "relSd", "relci95upper", "relci95lower") %in% colnames(resultSingle)))
  expect_equal(sum(resultSingle$relC),1)
})

test_that("single group results are the same as multi group results", {
  expect_equal(resultSingle$contrib, result$contrib[result$group=="A"], tolerance = 2)
  expect_equal(resultSingle$sdv, result$sdv[result$group=="A"], tolerance = 0.5)
  expect_equal(resultSingle$relDeviance, result$relDeviance[result$group=="A"])
})

