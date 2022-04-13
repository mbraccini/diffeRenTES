#New string append operator utilizzabile "new" %+% " string" o `%++%`("pippo" ,"pluto")
`%+%` <- function(a, b) paste0(a, b)
`%++%` <- function(a, b) paste(a, b, sep=" ")

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
    #NAMES
    namesAttrs <- c(1:length(attractors$attractors))
    names(syncAttractors) <- sapply("a", paste0,namesAttrs )
    return (syncAttractors)
}


#' Compute ATM
#'
#' \code{getATM} returns the ATM (Attractor Transistion Matrix) structure. 
#' The ATM computes the probability of a transition between the attractors of the Boolean network upon the introduction of noise in the form of a logic negation to each node of each state of each attractor,  checking in which attractor the dynamics relaxes.
#' The diagonal of the ATM accounts for attractor robustness, as diagonal values represent the probability of returning to the same attractor after a perturbation.
#'
#' @param net The Boolean network previously loaded with loadNetwork() of BoolNet package
#' @param syncAttractors Synchronous attractors of the Boolean network
#' @param MAX_STEPS_TO_FIND_ATTRACTORS Number of steps after that the dynamics after the perturbation gives up
#'
#' @return The output will be a named list containing the computed ATM structure, the number of the lost flips (i.e., the number of perturbations that have not reach another attractor within the provided MAX_STEPS_TO_FIND_ATTRACTORS), and lastly the attractors in two formats: the one returned by the BoolNet package (called decimal) and their binary translation (called binary).
#'
#' @import BoolNet
#' @import igraph
#' @examples
#' 
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net) 
#' getATM(net, attractors)
#'
#' @export
getATM <- function(net, syncAttractors, MAX_STEPS_TO_FIND_ATTRACTORS = 1000){
    initialAttractors <- syncAttractors
    numGenes <- length(syncAttractors$stateInfo$genes)
    numAttractors <- length(syncAttractors$attractors)
    
    attractors <- getMatricesAttractors(syncAttractors,numGenes)
    ATM <- matrix( rep( 0, len=numAttractors^2), nrow = numAttractors)
    #print(attractors[1][[1]][,1])
    
    lost <- 0
    for(i in c(1:numAttractors)){
        att <- attractors[i][[1]]
        for(j in 1:ncol(att)){ 
            for (gene in c(1:numGenes)){
                flipped <- att[,j]
                flipped[gene] <- xor(att[,j][gene], 1)
                attFoundIndex <- evolveUntilAttractor(net, attractors, flipped, MAX_STEPS_TO_FIND_ATTRACTORS)
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
    
    attrsDecimal <- initialAttractors$attractors
    names(attrsDecimal) <- names(attractors)
    attrs <- list("decimal" = attrsDecimal, "binary"=attractors)
    a <- list ("ATM" = ATM, "lostFLips" = lost, "attractors"= attrs)
    return (a)
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
evolveUntilAttractor <- function(net, attractors, state, MAX_STEPS_TO_FIND_ATTRACTORS){
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


#' Compute TES
#'
#' Creates a structure for constructing the TES as described in "A Dynamical Model of Genetic Networks for Cell Differentiation 
#'Villani M, Barbieri A, Serra R (2011) A Dynamical Model of Genetic Networks for Cell Differentiation. PLOS ONE 6(3): e17703. https://doi.org/10.1371/journal.pone.0017703"
#'
#' @param ATM ATM structure as returned from the \code{\link{getATM}} method
#'
#' @return The output will be a named list that contains the list of computed TESs, the noise thresholds at which they emerged and lastly the ATM structure.
#'
#' @examples
#' 
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net) 
#' ATM <- getATM(net, attractors)
#' getTESs(ATM) 
#'
#' @export
getTESs <- function(ATM){
    ATMstructure <- ATM
    ATM <- ATM$ATM
    thresholds <- unique(sort(c(ATM)))
    if (!(0 %in% thresholds)){
        print("There is no threshold value equal to 0 in the ATM, I will add it for the TESs computation")
        thresholds <- c(0, thresholds)
    }
    tes_i <- 1
    tes_list <- list()
    totalTES <- 0
    for (thrs in thresholds){
        
        ATM[ATM <= thrs] <- 0
        adjMtrx <- ATM
        adjMtrx[adjMtrx != 0] <- 1
       

        ATMgraph <- graph_from_adjacency_matrix(adjMtrx, mode="directed")
        scc <- clusters(ATMgraph, mode="strong")


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
        
        namesTESs <- (c(1:length(tes))) + totalTES
        names(tes) <- sapply("TES_", paste0,namesTESs)
        
        
        tes_list[[tes_i]] <- tes
        totalTES <- totalTES + length(tes_list[[tes_i]])
        
        tes_i <- tes_i + 1
    }
    
    namesLEVELTESs <- c(1:length(tes_list)) - 1
    names(tes_list) <- sapply("level_",  paste0, namesLEVELTESs)
    a <- list("TES"=tes_list, "thresholds"=thresholds, "ATM"=ATMstructure)
    return (a)
}



checkUpperlevels <- function(attrs, tesLvl, grandFatherLevel){
    if (grandFatherLevel < 1){
        return (-1)
    } else {
        for(grandFatherTESname in names(tesLvl[[grandFatherLevel]])){ #livello sopra per cercare i padri
            if (attrs %in% tesLvl[[grandFatherLevel]][[grandFatherTESname]]){
                return (grandFatherTESname)
            }
        }
        checkUpperlevels(attrs,  tesLvl, grandFatherLevel - 1)
    }
}


#' Save the differentiation tree using the DOT syntax and, if needed, produce the related image.
#'
#' \code{saveDifferentiationTreeToFile} saves the DOT representation of the derived differentiation tree into a file.
#'
#' @param TESs TES structure computed with \code{\link{getTESs}}
#' @param filename The filename of the .gv file
#' @param saveImage Logical parameter indicating whether \code{saveDotFileDifferentiationTree} have to produce also the image of the tree, in .svg format.
#'
#' @return None
#'
#' @import DOT
#'
#' @examples
#'
#' \dontrun{
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net) 
#' ATM <- getATM(net, attractors)
#' TESs <- getTESs(ATM) 
#' saveDifferentiationTreeToFile(TESs, "exampleTree")
#' }	
#'
#' @export
saveDifferentiationTreeToFile <- function(TESs, filename, saveImage=TRUE){
	DOTRep <- getDifferentiationTreeAsDOTString(TESs)
    sink(filename %+% ".gv")
    cat(DOTRep)
    sink()
    if (saveImage){
        dot(DOTRep, file = filename %+% ".svg")
    }
}


# Generates a DOT tree-like representation of the TESs
getDifferentiationTreeAsDOTString <- function(TESs){
    tesLvl <- TESs$TES
    namesPerLevel <- function(level){ names(level)}
    TESnames <- lapply(tesLvl, namesPerLevel)
    TESnames <- unlist(TESnames)
    gString <- "digraph diffTree {forcelabels=true;\n"
    for(lvl in c(1:length(tesLvl))){#per ogni livello
        if (lvl != 1) {#lvl 1=livello 0 dell'albero
            #Check if this level is equal to the previous one
            if (length(tesLvl[[lvl]]) == length(TESs$ATM$attractors$decimal)){ #if number of TES is the same of the number of attractors 
                if(length(tesLvl[[lvl]]) == length(tesLvl[[lvl - 1]])){
                    #same number of TES for this 2 adjacent levels
                    #isOne <- function(x) x == 1
                    #if (all(sapply(sapply(tesLvl[[lvl]],FUN=length),FUN=isOne)) && all(sapply(sapply(tesLvl[[lvl - 1]],FUN=length),FUN=isOne))){
                        next
                        #}
                }
            }
        }
        #ADDING NODES
        for(tesNAME in names(tesLvl[[lvl]])){#per ogni TES nel livello in esame
            attractorsOfThisTES <- "attrs:" %+% paste(tesLvl[[lvl]][[tesNAME]],collapse=",")
            gString <- gString %+% tesNAME %+%" [label = \"" %+% tesNAME %+% "\\n "%+% attractorsOfThisTES %+%"\"];\n"
            found <- FALSE
            if (lvl != 1) {#lvl 1=livello 0 dell'albero
                for(fatherTESname in names(tesLvl[[lvl -1]])){ #livello sopra per cercare i padri
                    if (tesLvl[[lvl]][[tesNAME]] %in% tesLvl[[lvl - 1]][[fatherTESname]]){
                        gString <- gString %+% fatherTESname %+% " -> " %+% tesNAME %+% "[label="%+% TESs$thresholds[[lvl]] %+%"];\n"
                        found <- TRUE
                    }                
                }
                if (!found) {
                    #check in upper levels
                    grandFatherTesNAME <- checkUpperlevels(tesLvl[[lvl]][[tesNAME]], tesLvl, lvl - 2)
                    if (grandFatherTesNAME != -1){
                        gString <- gString %+% grandFatherTesNAME %+% " -> " %+% tesNAME %+% "[style=dashed, color=grey];\n"
                    }
                }
            }
        }
        #ADDING RANKS
        sameRANK <- sapply(names(tesLvl[[lvl]]),  paste0, ";")
        sameRANK <- paste(sameRANK, collapse = "")
        gString <- gString %+% "{ rank=same;" %+% sameRANK %+% " }\n"
    }
    sString <- gString %+% "}"
    return (sString)
}
 

main <- function(){
    #fileBN <- "test/self_loop_bn_1_BoolNet.bn"
    #fileBN <- "rete.bn"
    #net <- loadNetwork(fileBN)
    net <- generateRandomNKNetwork(20, 2)
    #saveNetwork(net,"rete.bn")

    attractors <-  getAttractors(net) 

    ATM <- getATM(net, attractors)
    print(ATM)
    print(ATM$attractors$decimal$a2)
    print(ATM$attractors$binary$a2)


    print(ATM$ATM)

    TESs <- getTESs(ATM)

    print(TESs)

    #print(TESs$TES$"level_2"$"TES_4")
    strg <- getDifferentiationTreeAsDOTString(TESs)
    print(strg)
    saveDotFileDifferentiationTree(strg, "diffeTREE2")
}


