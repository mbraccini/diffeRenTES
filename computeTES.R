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



computeTESs <- function(ATM){
    thresholds <- unique(ATM)
    thresholds <- c(0, thresholds)
    tes_i <- 1
    tes_list <- list()
    
    
    for (thrs in thresholds){
        
        ATM[ATM <= thrs] <- 0
        adjMtrx <- ATM
        adjMtrx[adjMtrx != 0] <- 1
        #print(ATM)
        #print(adjMtrx)
        #adjMtrx <- matrix( rep( 0, len=25), nrow = 5)
        #adjMtrx[1,1] <- 1
        #adjMtrx[2,2] <- 1
        #adjMtrx[2,5] <- 1
        #adjMtrx[1,2] <- 1
        #adjMtrx[2,1] <- 1
        #adjMtrx[3,2] <- 1
        #adjMtrx[4,3] <- 1
        #rownames(adjMtrx) <- c("A1","A2","A3","A4","A5")
        #colnames(adjMtrx) <- c("A1","A2","A3","A4","A5")

        ATMgraph <- graph_from_adjacency_matrix(adjMtrx, mode="directed")
        #print(ATMgraph)
        scc <- clusters(ATMgraph, mode="strong")
        #plot(ATMgraph)


        scc_list <- list()
        i <- 1
        scc_unique_id <- unique(scc$membership)
        for (id in scc_unique_id){
            scc_list[[i]] <- names(which(scc$membership==id)) #prendo i nomi di tutti gli elm di quella scc
            i <- i+1
        }
    
        tes<- list()
        j <- 1
        for(scc in scc_list){
            isTES <- TRUE
            for(elm in scc){
                #print(scc)
                out_arcs <- names(which(adjMtrx[elm,] > 0))
                if (any(!(out_arcs %in% scc))){
                    isTES <- FALSE
                    break
                }
            }
            if (isTES){
                tes[[j]] <- scc                
                j <- j + 1
            }
        }
        namesTESs <- (c(1:length(tes)) - 1 ) + tes_i
        names(tes) <- sapply(namesTESs, paste0, "TES")
        
        
        tes_list[[tes_i]] <- tes
        tes_i <- tes_i + 1
    }
    
    namesLEVELTESs <- c(1:length(tes_list)) - 1
    names(tes_list) <- sapply(namesLEVELTESs, paste0, "level")
    a <- list("tes"=tes_list, "thresholds"=thresholds)
    return (a)
}

computeDiffTree <- function(TESs){
    tesLvl <- TESs$tes
    namesPerLevel <- function(level){ names(level)}
    TESnames <- lapply(tesLvl, namesPerLevel)
    TESnames <- unlist(TESnames)
    
    adjMatrixDiffTree <- matrix( rep( 0, len= length(TESnames)^2), nrow = length(TESnames))
    rownames(adjMatrixDiffTree) <- TESnames
    colnames(adjMatrixDiffTree) <- TESnames
    
    
    #adjMatrixDiffTree[1,4] <- 1
    #adjMatrixDiffTree[1,3] <- 1
    print(tesLvl)
    for(level in c(2:length(tesLvl))){
        print(paste0("level  ", level -1 ))
        #print(tesLvl[[level]])
        for(tesName in names(tesLvl[[level]])){
            #print(tesLvl[[level]][[tesName]])
           
        }
        
    }
    
    print(adjMatrixDiffTree)
    
    diffGraph <- graph_from_adjacency_matrix(adjMatrixDiffTree, mode="directed")
    #plot(diffGraph, layout = layout.reingold.tilford(diffGraph, root=names(tesLvl[[1]])))
    
    
    l <- layout_with_sugiyama(diffGraph)
    plot(l$extd_graph, 
         vertex.shape="circle", 
         vertex.label=as_ids(V(diffGraph)),
         edge.arrow.mode=2, 
         edge.arrow.width=2, 
         edge.arrow.size=0.1)
    
}

ATM <- computeATM("test/self_loop_bn_1_BoolNet.bn")
#print(ATM)
TESs <- computeTESs(ATM)
print(TESs$tes$"2level"$"4TES")
computeDiffTree(TESs)

