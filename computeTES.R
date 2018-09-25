library(BoolNet)
library(igraph)


MAX_STEPS_TO_FIND_ATTRACTORS <- 1000

#Function - from attractors and gene numbers retrieve a list of matrices with columns that are the states that compose the attractors of the net
getMatricesAttractors <- function(attractors, numGens){
    matrixAttractorList <- list()
    #print(attractors$attractors[[1]]$involvedStates) #A matrix containing the states that make up the attractor. Each column represents one state.  The entries are decimal numbers that internally represent the states of the genes
    for (i in c(1:length(attractors$attractors))){ #per il numero di attrattori
        matrixAttractorList[[i]] <- attractors$attractors[[i]]$involvedStates
    }
    attsNames <- attractors$stateInfo$genes  #ci appiccichiamo anche il nome dei geni
    intToBitVector <- function(attractorColumn, numBit){ 
        v <- as.numeric(intToBits(attractorColumn))[1:numBit]
        names(v) <- attsNames
        return (v)
    }
    funForAMatrixDescribingTheAttractor <- function(attractorsMatrix,numBits){ apply(attractorsMatrix,2,intToBitVector, numBit=numBits)} #2 significa applicalo sulle colonne
    syncAttractors <- lapply(matrixAttractorList,funForAMatrixDescribingTheAttractor,numBits=numGens) 
    #NOMI
    namesAttrs <- c(1:length(attractors$attractors))
    names(syncAttractors) <- sapply(namesAttrs, paste0, "A")
    return (syncAttractors)
}


computeATM <- function(fileBN){
    net <- loadNetwork(fileBN)
    numGenes <- length(net$genes)
    attractors <- getAttractors(net) 
    numAttractors <- length(attractors$attractors)
    attractors <- getMatricesAttractors(attractors,numGenes)
    
    ATM <- matrix( rep( 0, len=numAttractors^2), nrow = numAttractors)
    #print(attractors[1][[1]][,1])
    print(attractors)
    #print(attractors$"2att")
    
    lost <- 0
    for(i in c(1:numAttractors)){
        att <- attractors[i][[1]]
        for(j in 1:ncol(att)){ 
            for (gene in c(1:numGenes)){
                flipped <- att[,j]
                flipped[gene] <- xor(att[,j][gene], 1)
                attFoundIndex <- evolveUntilAttractor(net, attractors, flipped)
                if (attFoundIndex != -1){
                    ATM[i,attFoundIndex] <- ATM[i,attFoundIndex] + 1
                } else {
                    lost <- lost + 1
                }
            }
        }
    }
    
    rownames(ATM) <- names(attractors)
    colnames(ATM) <- names(attractors)
    
    ATM <- sweep(ATM, 1, rowSums(ATM), FUN="/") #normalization per row sums
    ATM <- round(ATM, 2)
    return (ATM)
}

matchAgainstKnownAttractors <- function(state, attractors){
    numAttractors <- length(attractors)
    for(i in c(1:numAttractors)){
        att <- attractors[i][[1]]
        for(j in 1:ncol(att)){ 
            if (all(state==att[,j])){
                return (i)
            }
        }
    }
    return (-1)
}

#Returns index of attractors found or -1 if too much time needed
evolveUntilAttractor <- function(net, attractors, state){
    for (i in c(1:MAX_STEPS_TO_FIND_ATTRACTORS)){
        idx <- matchAgainstKnownAttractors(state,attractors)
        if(idx != -1){
            return (idx)
        } else {
            state <- stateTransition(net, state, type=c("synchronous"))
        }
    }
    return (-1)
}


ATM <- computeATM("test/self_loop_bn_1_BoolNet.bn")
ATM[ATM < 0.8] <- 0
adjMtrx <- ATM
adjMtrx[adjMtrx != 0] <- 1
print(adjMtrx)
print(ATM)

ATMgraph <- graph_from_adjacency_matrix(adjMtrx, mode="directed")
print(ATMgraph)
scc <- clusters(ATMgraph, mode="strong")

print(scc)
