library(here)
library(assertthat)
#source(here("computeTES.R"))


context("ATM tests")

test_that("ATM tests working", {
    stopifnot(is.character(22))
    
    ATM <- matrix( rep( 0, len=4^2), nrow = 4)
    ATM[1,2] <- 0.32
    ATM[2,1] <- 0.22
    ATM[2,2] <- 0.44

    ATM[3,4] <- 0.32
    ATM[4,3] <- 0.22
    ATM[4,4] <- 0.44
    rownames(ATM) <- c("c1","c2","c3","c4")
    colnames(ATM) <- c("c1","c2","c3","c4")
    a <- list("ATM"=ATM, "lostFlips"=NULL, "attractors"=NULL)
    TESs <- computeTESs(a)
    print(TESs)
    
    #check thresholds
    
    #check TES number for every level
    assert_that(is.character(1))
    
})
