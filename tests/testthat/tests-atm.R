#' @importFrom BoolNet getAttractors
context("TES tests")

test_that("TES tests working", {
  ATM <- matrix(rep(0, len = 4^2), nrow = 4)
  ATM[1, 2] <- 0.32
  ATM[2, 1] <- 0.22
  ATM[2, 2] <- 0.44

  ATM[3, 4] <- 0.32
  ATM[4, 3] <- 0.22
  ATM[4, 4] <- 0.44
  rownames(ATM) <- c("a1", "a2", "a3", "a4")
  colnames(ATM) <- c("a1", "a2", "a3", "a4")

  a <- list("ATM" = ATM, "lostFlips" = NULL, "attractors" = NULL)
  TESs <- getTESs(a)

  expect_equal(TESs$thresholds, c(0, 0.22, 0.32, 0.44))

  first <- TESs$TES[[1]]
  expect_equal(length(first), 2) # with first threshold we must have 2 TESs
  for (t in first) {
    expect_equal(length(t), 2) # each first threshold TES must have 2 attractors
  }
  second <- TESs$TES[[2]]
  expect_equal(length(second), 2) # with second threshold we must have 2 TESs
  for (t in second) {
    expect_equal(length(t), 1) # each second threshold TES must have 1 attractor
  }
  third <- TESs$TES[[3]]
  expect_equal(length(third), 4) # with third threshold we must have 4 TESs
  for (t in third) {
    expect_equal(length(t), 1) # each third threshold TES must have 1 attractor
  }
  fourth <- TESs$TES[[4]]
  expect_equal(length(fourth), 4) # with fourth threshold we must have 4 TESs
  for (t in fourth) {
    expect_equal(length(t), 1) # each fourth threshold TES must have 1 attractor
  }
})


test_that("test ATM rete self_loop_bn_1_BoolNet.bn", {
  n <- BoolNet::loadNetwork("self_loop_bn_1_BoolNet.bn")
  res <- diffeRenTES::getATM(net = n, synchronous_attractors = BoolNet::getAttractors(n))
  ATM <- res$ATM
  expect_equal(ATM[1, 1], 0.33)
  expect_equal(ATM[1, 2], 0.67)
  expect_equal(ATM[2, 1], 0.33)
  expect_equal(ATM[2, 2], 0.67)
})


test_that("test ATM rete self_loop_bn_1_BoolNet.bn", {
  n <- BoolNet::loadNetwork("attractor_landscape_paper.bn")
  res <- diffeRenTES::getATM(net = n, synchronous_attractors = BoolNet::getAttractors(n))
  ATM <- res$ATM
  expect_equal(ATM[1, 1], 0)
  expect_equal(ATM[1, 2], 0.75)
  expect_equal(ATM[1, 3], 0.25)
  expect_equal(ATM[1, 4], 0)

  expect_equal(ATM[2, 1], 0.25)
  expect_equal(ATM[2, 2], 0)
  expect_equal(ATM[2, 3], 0.25)
  expect_equal(ATM[2, 4], 0.5)

  expect_equal(ATM[3, 1], 0)
  expect_equal(ATM[3, 2], 0)
  expect_equal(ATM[3, 3], 0.25)
  expect_equal(ATM[3, 4], 0.75)

  expect_equal(ATM[4, 1], 0)
  expect_equal(ATM[4, 2], 0.5)
  expect_equal(ATM[4, 3], 0.38)
  expect_equal(ATM[4, 4], 0.12)
})
