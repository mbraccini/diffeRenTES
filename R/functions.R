`%+%` <- function(a, b) paste0(a, b)
`%++%` <- function(a, b) paste(a, b, sep = " ")

get_matrices_attractors <- function(attractors, numGens) {
  matrix_attractor_list <- list()
  for (i in c(1:length(attractors$attractors))) {
    matrix_attractor_list[[i]] <- attractors$attractors[[i]]$involvedStates
  }
  atts_names <- attractors$stateInfo$genes
  int_to_bit_vector <- function(attractor_column, num_bit) {
    v <- as.numeric(intToBits(attractor_column))[1:num_bit]
    names(v) <- atts_names
    return(v)
  }
  fun_for_a_matrix_describing_the_attractor <- function(attractors_matrix, num_bits) {
    apply(attractors_matrix, 2, int_to_bit_vector, num_bit = num_bits)
  }
  synchronous_attractors <- lapply(matrix_attractor_list, fun_for_a_matrix_describing_the_attractor, num_bits = numGens)
  names_attrs <- c(1:length(attractors$attractors))
  names(synchronous_attractors) <- sapply("a", paste0, names_attrs)
  return(synchronous_attractors)
}


#' Compute ATM
#'
#' \code{getATM} returns the ATM (Attractor Transistion Matrix) structure.
#' The ATM computes the probability of a transition between the attractors of the Boolean network upon the introduction of noise in the form of a logic negation to each node of each state of each attractor,  checking in which attractor the dynamics relaxes.
#' The diagonal of the ATM accounts for attractor robustness, as diagonal values represent the probability of returning to the same attractor after a perturbation.
#'
#' @param net The Boolean network previously loaded with loadNetwork() of BoolNet package
#' @param synchronous_attractors Synchronous attractors of the Boolean network
#' @param MAX_STEPS_TO_FIND_ATTRACTORS Number of steps after that the dynamics after the perturbation gives up
#'
#' @return The output will be a named list containing the computed ATM structure, the number of the lost flips (i.e., the number of perturbations that have not reach another attractor within the provided MAX_STEPS_TO_FIND_ATTRACTORS), and lastly the attractors in two formats: the one returned by the BoolNet package (called decimal) and their binary translation (called binary).
#'
#' @importFrom BoolNet stateTransition
#' @examples
#'
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net)
#' getATM(net, attractors)
#'
#' @export
getATM <- function(net, synchronous_attractors, MAX_STEPS_TO_FIND_ATTRACTORS = 1000) {
  initial_attractors <- synchronous_attractors
  num_genes <- length(synchronous_attractors$stateInfo$genes)
  num_attractors <- length(synchronous_attractors$attractors)

  attractors <- get_matrices_attractors(synchronous_attractors, num_genes)
  ATM <- matrix(rep(0, len = num_attractors^2), nrow = num_attractors)

  lost <- 0
  for (i in c(1:num_attractors)) {
    att <- attractors[i][[1]]
    for (j in 1:ncol(att)) {
      for (gene in c(1:num_genes)) {
        flipped <- att[, j]
        flipped[gene] <- xor(att[, j][gene], 1)
        att_found_index <- evolve_until_attractor(net, attractors, flipped, MAX_STEPS_TO_FIND_ATTRACTORS)
        if (att_found_index != -1) {
          ATM[i, att_found_index] <- ATM[i, att_found_index] + 1
        } else {
          lost <- lost + 1
        }
      }
    }
  }

  rownames(ATM) <- names(attractors)
  colnames(ATM) <- names(attractors)

  ATM <- sweep(ATM, 1, rowSums(ATM), FUN = "/")
  ATM <- round(ATM, 2)

  attrs_decimal <- initial_attractors$attractors
  names(attrs_decimal) <- names(attractors)
  attrs <- list("decimal" = attrs_decimal, "binary" = attractors)
  a <- list("ATM" = ATM, "lostFLips" = lost, "attractors" = attrs)
  return(a)
}

match_against_known_attractors <- function(state, attractors) {
  num_attractors <- length(attractors)
  for (i in c(1:num_attractors)) {
    att <- attractors[i][[1]]
    for (j in 1:ncol(att)) {
      if (all(state == att[, j])) {
        return(i)
      }
    }
  }
  return(-1)
}

evolve_until_attractor <- function(net, attractors, state, MAX_STEPS_TO_FIND_ATTRACTORS) {
  for (i in c(1:MAX_STEPS_TO_FIND_ATTRACTORS)) {
    idx <- match_against_known_attractors(state, attractors)
    if (idx != -1) {
      return(idx)
    } else {
      state <- stateTransition(net, state, type = c("synchronous"))
    }
  }
  return(-1)
}

#' Compute TES
#'
#' Creates a structure for constructing the TES as described in "A Dynamical Model of Genetic Networks for Cell Differentiation
#' Villani M, Barbieri A, Serra R (2011) A Dynamical Model of Genetic Networks for Cell Differentiation. PLOS ONE 6(3): e17703. https://doi.org/10.1371/journal.pone.0017703"
#'
#' @param ATM ATM structure as returned from the \code{\link{getATM}} method
#'
#' @return The output will be a named list that contains the list of computed TESs, the noise thresholds at which they emerged and lastly the ATM structure.
#'
#' @importFrom igraph graph_from_adjacency_matrix clusters
#'
#' @examples
#'
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net)
#' ATM <- getATM(net, attractors)
#' getTESs(ATM)
#'
#' @export
getTESs <- function(ATM) {
  ATM_structure <- ATM
  ATM <- ATM$ATM
  thresholds <- unique(sort(c(ATM)))
  if (!(0 %in% thresholds)) {
    print("There is no threshold value equal to 0 in the ATM, I will add it for the TESs computation")
    thresholds <- c(0, thresholds)
  }
  tes_i <- 1
  tes_list <- list()
  total_tes <- 0
  for (thrs in thresholds) {
    ATM[ATM <= thrs] <- 0
    adj_mtrx <- ATM
    adj_mtrx[adj_mtrx != 0] <- 1


    ATM_graph <- igraph::graph_from_adjacency_matrix(adj_mtrx, mode = "directed")
    scc <- igraph::clusters(ATM_graph, mode = "strong")


    scc_list <- list()
    i <- 1
    scc_unique_id <- unique(scc$membership)
    for (id in scc_unique_id) {
      scc_list[[i]] <- names(which(scc$membership == id))
      i <- i + 1
    }

    tes <- list()
    j <- 1
    for (scc in scc_list) {
      is_TES <- TRUE
      for (elm in scc) {
        out_arcs <- names(which(adj_mtrx[elm, ] > 0))
        if (any(!(out_arcs %in% scc))) {
          is_TES <- FALSE
          break
        }
      }
      if (is_TES) {
        tes[[j]] <- scc
        j <- j + 1
      }
    }

    names_TESs <- (c(1:length(tes))) + total_tes
    names(tes) <- sapply("TES_", paste0, names_TESs)


    tes_list[[tes_i]] <- tes
    total_tes <- total_tes + length(tes_list[[tes_i]])

    tes_i <- tes_i + 1
  }

  names_level_TESs <- c(1:length(tes_list)) - 1
  names(tes_list) <- sapply("level_", paste0, names_level_TESs)
  a <- list("TES" = tes_list, "thresholds" = thresholds, "ATM" = ATM_structure)
  return(a)
}



check_upperlevels <- function(attrs, tes_lvl, grand_father_level) {
  if (grand_father_level < 1) {
    return(-1)
  } else {
    for (grand_father_TES_name in names(tes_lvl[[grand_father_level]])) {
      if (attrs %in% tes_lvl[[grand_father_level]][[grand_father_TES_name]]) {
        return(grand_father_TES_name)
      }
    }
    check_upperlevels(attrs, tes_lvl, grand_father_level - 1)
  }
}


#' Save the differentiation tree using the DOT syntax and, if needed, produce the related image.
#'
#' \code{saveDifferentiationTreeToFile} saves the DOT representation of the derived differentiation tree into a file.
#'
#' @param TESs TES structure computed with \code{\link{getTESs}}
#' @param filename The filename of the .gv file
#' @param save_image Logical parameter indicating whether \code{saveDotFileDifferentiationTree} have to produce also the image of the tree, in .svg format.
#'
#' @return None
#'
#' @importFrom DOT dot
#'
#' @examples
#' \dontrun{
#' net <- BoolNet::generateRandomNKNetwork(10, 2)
#' attractors <- BoolNet::getAttractors(net)
#' ATM <- getATM(net, attractors)
#' TESs <- getTESs(ATM)
#' saveDifferentiationTreeToFile(TESs, "exampleTree")
#' }
#'
#' @export
saveDifferentiationTreeToFile <- function(TESs, filename, save_image = TRUE) {
  dot_rep <- get_tree_as_dot_string(TESs)
  sink(filename %+% ".gv")
  cat(dot_rep)
  sink()
  if (save_image) {
    dot(dot_rep, file = filename %+% ".svg")
  }
}

get_tree_as_dot_string <- function(TESs) {
  tes_lvl <- TESs$TES
  names_per_level <- function(level) {
    names(level)
  }
  TES_names <- lapply(tes_lvl, names_per_level)
  TES_names <- unlist(TES_names)
  g_string <- "digraph diffTree {forcelabels=true;\n"
  for (lvl in c(1:length(tes_lvl))) {
    if (lvl != 1) {
      if (length(tes_lvl[[lvl]]) == length(TESs$ATM$attractors$decimal)) {
        if (length(tes_lvl[[lvl]]) == length(tes_lvl[[lvl - 1]])) {
          next
        }
      }
    }

    for (tesNAME in names(tes_lvl[[lvl]])) {
      attractors_of_this_TES <- "attrs:" %+% paste(tes_lvl[[lvl]][[tesNAME]], collapse = ",")
      g_string <- g_string %+% tesNAME %+% " [label = \"" %+% tesNAME %+% "\\n " %+% attractors_of_this_TES %+% "\"];\n"
      found <- FALSE
      if (lvl != 1) {
        for (fatherTESname in names(tes_lvl[[lvl - 1]])) {
          if (tes_lvl[[lvl]][[tesNAME]] %in% tes_lvl[[lvl - 1]][[fatherTESname]]) {
            g_string <- g_string %+% fatherTESname %+% " -> " %+% tesNAME %+% "[label=" %+% TESs$thresholds[[lvl]] %+% "];\n"
            found <- TRUE
          }
        }
        if (!found) {
          grand_father_TES_name <- check_upperlevels(tes_lvl[[lvl]][[tesNAME]], tes_lvl, lvl - 2)
          if (grand_father_TES_name != -1) {
            g_string <- g_string %+% grand_father_TES_name %+% " -> " %+% tesNAME %+% "[style=dashed, color=grey];\n"
          }
        }
      }
    }

    same_rank <- sapply(names(tes_lvl[[lvl]]), paste0, ";")
    same_rank <- paste(same_rank, collapse = "")
    g_string <- g_string %+% "{ rank=same;" %+% same_rank %+% " }\n"
  }
  s_string <- g_string %+% "}"
  return(s_string)
}
