#!/usr/bin/env Rscript

library(CellNOptR)
library(here)
source(here::here("tools", "functions.R"))

#' Perturb a Boolean-logic model by random local edits
#' 
#' @param model  A CellNOptR model list (interMat, notMat, namesSpecies, reacID)
#' @param change_percent  A number in [0,1].  Fraction of nodes to perturb.
#' @return  A modified copy of `model`

perturbModel <- function(model, cnolist, change_percent=0.1) {
  # Helper to add a new reaction column
  add_reaction <- function(mod, new_id, src, tgt, is_not=FALSE) {
    # skip if already present
    if (new_id %in% mod$reacID) return(mod)
    # add to reacID
    mod$reacID <- c(mod$reacID, new_id)
    # extend interMat and notMat
    nsp <- length(mod$namesSpecies)
    # add zero column
    mod$interMat <- cbind(mod$interMat, rep(0, nsp))
    mod$notMat   <- cbind(mod$notMat,   rep(0, nsp))
    col_idx <- ncol(mod$interMat)
    # source → -1, target → +1
    mod$interMat[src, col_idx] <- -1
    mod$interMat[tgt, col_idx] <- +1
    # mark NOT if needed
    if (is_not) mod$notMat[src, col_idx] <- 1

    # **re‑assign column names**
    colnames(mod$interMat) <- mod$reacID
    colnames(mod$notMat)   <- mod$reacID
    return(mod)
  }
  
  # Helper to remove a reaction column by index
  remove_reaction <- function(mod, col_idx) {
    mod$reacID   <- mod$reacID[-col_idx]
    mod$interMat <- mod$interMat[ , -col_idx]
    mod$notMat   <- mod$notMat[ , -col_idx]

    # **re‑assign column names**
    colnames(mod$interMat) <- mod$reacID
    colnames(mod$notMat)   <- mod$reacID
    return(mod)
  }
  
  # Parse a reaction string "(!)node1=node2" into components
  parse_reac <- function(rid) {
    neg <- startsWith(rid, "!")
    rid2 <- sub("^!", "", rid)
    parts <- strsplit(rid2, "=", fixed=TRUE)[[1]]
    list(src=parts[1], tgt=parts[2], is_not=neg)
  }

  check_and_reverse_stimuli_edges <- function(src, tgt, cnolist) {
    stimuli <- cnolist$namesStimuli
    message(sprintf("Checking edge: %s → %s in %s", src, tgt, stimuli))
    # If target is a stimulus and source is not, reverse the edge
    if (tgt %in% stimuli && !(src %in% stimuli)) {
      # Reverse edge: stimulus should be source, not target
      message(sprintf("Reversing edge: %s → %s", tgt, src))
      return(list(src = tgt, tgt = src, reversed = TRUE))
    }
    # If both are stimuli, suggest connecting to a signal instead
    if (src %in% stimuli && tgt %in% stimuli) {
      signal <- sample(cnolist$namesSignals, 1)
      message(sprintf("Connecting stimuli %s and %s to signal %s", src, tgt, signal))
      return(list(src = src, tgt = signal, replaced = TRUE))
    }
    # Otherwise, keep edge as is
    return(list(src = src, tgt = tgt, unchanged = TRUE))
  }

  n_nodes <- length(model$namesSpecies)
  n_perturb <- ceiling(change_percent * n_nodes)
  sel_nodes <- sample(n_nodes, n_perturb)
  
  message(sprintf("Perturbing %d nodes out of %d (%.0f%%)", n_perturb, n_nodes, change_percent * 100))
  for (i in sel_nodes) {
    node1 <- model$namesSpecies[i]
    x <- runif(1)
    
    ## 1) Operator insertion (≈1/3 of the time)
    if (x < 1/3) {
      # pick random other node
      message(sprintf("Inserting operator for node %s", node1))
      j <- sample(setdiff(seq_len(n_nodes), i), 1)
      node2 <- model$namesSpecies[j]
      # choose operator type
      op <- sample(c("NOT","AND","OR"),1)
      if (op == "NOT") {
        message(sprintf("Working on case 1, subcase A"))
        rel <- check_and_reverse_stimuli_edges(node1, node2, cnolist)
        new_id <- sprintf("!%s=%s", rel$src, rel$tgt)
        i <- which(model$namesSpecies == rel$src)
        j <- which(model$namesSpecies == rel$tgt)
        model <- add_reaction(model, new_id, src=i, tgt=j, is_not=TRUE)
      } else {
        message(sprintf("Working on case 1, subcase B"))
        rel <- check_and_reverse_stimuli_edges(node1, node2, cnolist)        
        new_id <- sprintf("%s=%s", rel$src, rel$tgt)
        i <- which(model$namesSpecies == rel$src)
        j <- which(model$namesSpecies == rel$tgt)
        model <- add_reaction(model, new_id, src=i, tgt=j, is_not=FALSE)
      }
      
    ## 2) Subtree deletion (≈1/3)
    } else if (x < 2/3) {
      message(sprintf("Deleting subtree for node %s", node1))
      # find all reactions involving node1
      cols <- which(grepl(node1, model$reacID, fixed=TRUE))
      message(sprintf("Found %d reactions [%s] for node %s", length(cols), paste(model$reacID[cols], collapse=", "), node1))
      if (length(cols) > 1) {
          col <- sample(cols, 1)  # pick one at random
          message(sprintf("Deleting reaction %s", model$reacID[col]))
          rid <- model$reacID[col]
          model <- remove_reaction(model, col)
          # if (startsWith(rid, "!")) {
          #     model$notMat[which(model$namesSpecies == pr$src), col] <- 0
          # } else {
          #     pr <- parse_reac(rid)
          #     model$interMat[which(model$namesSpecies == pr$src), col] <- 0
          #     model$interMat[which(model$namesSpecies == pr$tgt), col] <- 0
          # }
      }
      
    ## 3) Operator exchange (≈1/3)
    } else {
      message(sprintf("Exchanging operator for node %s", node1))

      # pick inner coin
      xi <- runif(1)
      pos_cols <- which(!startsWith(model$reacID, "!") & startsWith(model$reacID, node1))
      neg_cols <- which(startsWith(model$reacID, "!") & startsWith(sub("^!", "", model$reacID), node1))
      if (length(pos_cols) >= 1 && length(neg_cols) >= 1) {

        # Case 1: Both positive and negative relations exist
        c_pos <- sample(pos_cols,1)
        c_neg <- sample(neg_cols,1)
        rid_pos <- model$reacID[c_pos]
        rid_neg <- model$reacID[c_neg]          
        # remove both
        model <- remove_reaction(model, c_neg)
        model <- remove_reaction(model, c_pos)
        # rel1 ← "!" + rel2
        new_rel1 <- paste0("!", rid_pos)
        # rel2 ← "rel1"  (i.e. the positive of rid_neg)
        new_rel2 <- sub("^!", "", rid_neg)
        # insert them
        pr1 <- parse_reac(rid_pos)  # first node info
        pr2 <- parse_reac(rid_neg)  # second node info
        # new_rel1 on first node
        src1 <- which(model$namesSpecies==pr1$src)
        tgt1 <- which(model$namesSpecies==pr1$tgt)
        model <- add_reaction(model, new_rel1, src=src1, tgt=tgt1, is_not=TRUE)
        # new_rel2 on first↔second of original neg
        src2 <- which(model$namesSpecies==pr2$src)
        tgt2 <- which(model$namesSpecies==pr2$tgt)
        model <- add_reaction(model, new_rel2, src=src2, tgt=tgt2, is_not=FALSE)
      
      } else if (length(pos_cols) >= 1 && length(neg_cols) == 0) {
        # Case 2: Only positive relations exist
        c_pos <- sample(pos_cols, 1)
        rid_pos <- model$reacID[c_pos]
        model <- remove_reaction(model, c_pos)
        new_rel <- paste0("!", rid_pos)
        pr <- parse_reac(rid_pos)
        src <- which(model$namesSpecies == pr$src)
        tgt <- which(model$namesSpecies == pr$tgt)
        model <- add_reaction(model, new_rel, src=src, tgt=tgt, is_not=TRUE)
        
      } else if (length(pos_cols) == 0 && length(neg_cols) >= 1) {
        # Case 3: Only negative relations exist
        c_neg <- sample(neg_cols, 1)
        rid_neg <- model$reacID[c_neg]
        model <- remove_reaction(model, c_neg)
        new_rel <- sub("^!", "", rid_neg)
        pr <- parse_reac(rid_neg)
        src <- which(model$namesSpecies == pr$src)
        tgt <- which(model$namesSpecies == pr$tgt)
        model <- add_reaction(model, new_rel, src=src, tgt=tgt, is_not=FALSE)
      }
    }
  }
  return(model)
}


#------------------------------------------------------------------------------#
#   Function: runPerturbPipeline                                                #
#------------------------------------------------------------------------------#
#' @param orig_model   A CellNOptR model (list with interMat, notMat, namesSpecies, reacID).
#'                     Defaults to ToyModel.
#' @param CNOlist      A CNOlist object for SIF→BoolNet conversion. Defaults to CNOlistToy.
#' @param change_pct   Fraction in [0,1] of nodes to perturb. Default: 0.2.
#' @param seed         Integer RNG seed for reproducibility. Default: 42.
#' @param out_dir_sif  Directory to save the perturbed SIF & RData. Default: "output/caspo".
#' @param out_dir_boolnet Directory to save the BoolNet file. Default: "output/boolnet".
#' @return              A list with elements `mod_model`, `sif_fname`, `rdata_fname`, `boolnet_fname`.
#------------------------------------------------------------------------------#
runPerturbPipeline <- function(orig_model     = { data("ToyModel", package="CellNOptR"); ToyModel },
                               cnolist        = { data("CNOlistToy", package="CellNOptR"); CNOlistToy },
                               change_pct     = 0.9,
                               seed           = 42,
                               out_dir_sif    = "output/caspo",
                               out_dir_boolnet= "output/boolnet") {
  # 1) Perturb
  set.seed(seed)
  cat(">> DEBUG: number of nodes available = ", length(orig_model$namesSpecies), "\n")
  cat(">> DEBUG: change_pct = ", change_pct, "\n")
  wanted_size <- floor(change_pct * length(orig_model$namesSpecies))
  cat(">> DEBUG: computed size = ", wanted_size, "\n")
  
  mod_model <- perturbModel(orig_model, cnolist, change_pct)
  message(sprintf("Reactions before: %d", length(orig_model$reacID)))
  message(sprintf("Reactions after : %d\n",  length(mod_model$reacID)))
  
  # 2) Prepare output dirs
  dir.create(out_dir_sif,      showWarnings=FALSE, recursive=TRUE)
  dir.create(out_dir_boolnet,  showWarnings=FALSE, recursive=TRUE)
  
  # Filenames
  pct_lbl      <- sprintf("%.0f", change_pct*100)
  sif_fname    <- file.path(out_dir_sif,     sprintf("ModifiedToyModel_%s.sif",  pct_lbl))
  rdata_fname  <- file.path(out_dir_sif,     sprintf("ModifiedToyModel_%s.RData",pct_lbl))
  boolnet_fname<- file.path(out_dir_boolnet, sprintf("ModifiedToyModel_%s.bnet",  pct_lbl))
  
  # Print the result of ifelse for each reaction
  for (i in seq_along(mod_model$reacID)) {
    res <- ifelse(any(mod_model$notMat[, i] == 1), -1, 1)
    cat(sprintf("Reaction %s: %d\n", mod_model$reacID[i], res))
  }
  # 3) Save RData & SIF
  message(mod_model)
  save(mod_model, file=rdata_fname)
  writeSIF(mod_model, file=sif_fname, overwrite=TRUE)
  message("Wrote:\n - RData → ", rdata_fname, "\n - SIF   → ", sif_fname, "\n")
  
  # 4) Convert to BoolNet format
  SIFToBoolNet(sifFile     = sif_fname,
               boolnetFile = boolnet_fname,
               CNOlist     = cnolist,
               model       = orig_model,
               fixInputs   = FALSE,
               preprocess  = TRUE,
               ignoreAnds  = TRUE)
  message("BoolNet file written to: ", boolnet_fname)
  
  # Return paths & perturbed model
  invisible(list(mod_model     = mod_model,
                 sif_fname     = sif_fname,
                 rdata_fname   = rdata_fname,
                 boolnet_fname = boolnet_fname))
}

#--- Load required libraries -------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
})

#--- Define command-line options ---------------------------------------------
option_list <- list(
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="Path to original model RData or RDS file", metavar="FILE"),
  make_option(c("-c", "--cno"), type="character", default=NULL,
              help="Path to CNOlist object (RDS)", metavar="FILE"),
  make_option(c("-p", "--changePCT"), type="double", default=0.9,
              help="Change percentage [default %default]", metavar="DOUBLE"),
  make_option(c("-s", "--seed"), type="integer", default=44,
              help="Random seed [optional]", metavar="INT"),
  make_option(c("--outSIF"), type="character", default="outputs/sif",
              help="Output directory for .sif files [default %default]", metavar="DIR"),
  make_option(c("--outBoolNet"), type="character", default="outputs/boolnet",
              help="Output directory for BoolNet files [default %default]", metavar="DIR")
)

parser <- OptionParser(option_list=option_list,
                       description = "Run perturbation pipeline for a Boolean model")
opt <- parse_args(parser)

#--- Validate inputs ----------------------------------------------------------
# if (is.null(opt$model)) {
#   stop("Error: --model must be provided (path to your original model file).")
# }
# if (is.null(opt$cno)) {
#   stop("Error: --cno must be provided (path to your CNOlist file).")
# }

#--- Load input objects ------------------------------------------------------
# Assumes these files contain the needed R objects:
#   - model loaded into variable 'myCustomModel'
#   - CNOlist loaded into 'myCNOlist'
# message("Loading model from: ", opt$model)
# load(opt$model)        # e.g. loads `myCustomModel`
# message("Loading CNOlist from: ", opt$cno)
# myCNOlist <- readRDS(opt$cno)

#--- (Optional) set seed -----------------------------------------------------
if (!is.null(opt$seed)) {
  set.seed(opt$seed)
  message("Random seed set to: ", opt$seed)
}

#--- Run the pipeline --------------------------------------------------------
if (is.null(opt$model)) {  
  results <- runPerturbPipeline(
    change_pct      = opt$changePCT,
  )
} else {    
  results <- runPerturbPipeline(
    orig_model      = myCustomModel,
    CNOlist         = myCNOlist,
    change_pct      = opt$changePCT,
  )
}

