library(CellNOptR)
library(here)
source(here::here("tools", "functions.R"))

#' Perturb a Boolean-logic model by random local edits
#' 
#' @param model  A CellNOptR model list (interMat, notMat, namesSpecies, reacID)
#' @param change_percent  A number in [0,1].  Fraction of nodes to perturb.
#' @return  A modified copy of `model`

perturbModel <- function(model, change_percent=0.1) {
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

  n_nodes <- length(model$namesSpecies)
  n_perturb <- ceiling(change_percent * n_nodes)
  sel_nodes <- sample(n_nodes, n_perturb)
  cat("Perturbing %d nodes out of %d (%.0f%%)\n", n_perturb, n_nodes, change_percent * 100)
  for (i in sel_nodes) {
    node1 <- model$namesSpecies[i]
    x <- runif(1)
    
    ## 1) Operator insertion (≈1/3 of the time)
    if (x < 1/3) {
      # pick random other node
      cat("Inserting operator for node %s\n", node1)
      j <- sample(setdiff(seq_len(n_nodes), i), 1)
      node2 <- model$namesSpecies[j]
      # choose operator type
      op <- sample(c("NOT","AND","OR"),1)
      if (op == "NOT") {
        new_id <- sprintf("!%s=%s", node1, node2)
        model <- add_reaction(model, new_id, src=i, tgt=j, is_not=TRUE)
      } else {
        new_id <- sprintf("%s=%s", node1, node2)
        model <- add_reaction(model, new_id, src=i, tgt=j, is_not=FALSE)
      }
      
    ## 2) Subtree deletion (≈1/3)
    } else if (x < 2/3) {
      cat("Deleting subtree for node %s\n", node1)
      # find all reactions involving node1
      cols <- which(grepl(node1, model$reacID, fixed=TRUE))

      if (length(cols)>0) {
        col <- sample(cols, 1)
        rid <- model$reacID[col]
        pr  <- parse_reac(rid)
        # if it's a NOT-edge, zero out notMat; else zero interMat
        if (pr$is_not && model$notMat[which(model$namesSpecies==pr$src),col]==1) {
          # do nothing
        } else {
          # model <- remove_reaction(model, col)
          # If I delete a reaction, there is chance that the model is inconsistent 
          # after I save the SIF model. 
          # So I'd rather do nothing here.
        }
      }
      
    ## 3) Operator exchange (≈1/3)
    } else {
      cat("Exchanging operator for node %s\n", node1)

      # pick inner coin
      xi <- runif(1)
      
      # --- subcase A: xi < 0.5 ---------------------------------------------
      if (xi < 0.5) {
        # pick one positive relation (no "!") and one negative ("!")
        pos_cols <- which(!startsWith(model$reacID, "!"))
        neg_cols <- which(startsWith(model$reacID, "!"))
        if (length(pos_cols)>0 && length(neg_cols)>0) {
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
        }
        
      # --- subcase B: xi >= 0.5 --------------------------------------------
      } else {
        # pick a normal relation (no "!")
        cols <- which(!startsWith(model$reacID,"!"))
        if (length(cols)>0) {
          col <- sample(cols,1)
          rid <- model$reacID[col]
          pr  <- parse_reac(rid)
          # remove it
          model <- remove_reaction(model, col)
          # decide swap direction based on src/tgt
          if (startsWith(pr$src, "1")) {
            # src was of form "1", so we treat via NOT-mat
            new_r <- sprintf("!%s=%s", pr$tgt, pr$src)
            j <- which(model$namesSpecies==pr$tgt)
            k <- which(model$namesSpecies==pr$src)
            model <- add_reaction(model, new_r, src=j, tgt=k, is_not=TRUE)
          } else {
            # simple edge swap
            new_r <- sprintf("%s=%s", pr$tgt, pr$src)
            j <- which(model$namesSpecies==pr$tgt)
            k <- which(model$namesSpecies==pr$src)
            model <- add_reaction(model, new_r, src=j, tgt=k, is_not=FALSE)
          }
        }
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
                               CNOlist        = { data("CNOlistToy", package="CellNOptR"); CNOlistToy },
                               change_pct     = 0.2,
                               seed           = 42,
                               out_dir_sif    = "output/caspo",
                               out_dir_boolnet= "output/boolnet") {
  # 1) Perturb
  set.seed(seed)
  mod_model <- perturbModel(orig_model, change_pct)
  message(sprintf("Reactions before: %d", length(orig_model$reacID)))
  message(sprintf("Reactions after : %d\n",  length(mod_model$reacID)))
  
  # 2) Prepare output dirs
  dir.create(out_dir_sif,      showWarnings=FALSE, recursive=TRUE)
  dir.create(out_dir_boolnet,  showWarnings=FALSE, recursive=TRUE)
  
  # Filenames
  pct_lbl      <- sprintf("%.0f%%", change_pct*100)
  sif_fname    <- file.path(out_dir_sif,     sprintf("ModifiedToyModel_%s.sif",  pct_lbl))
  rdata_fname  <- file.path(out_dir_sif,     sprintf("ModifiedToyModel_%s.RData",pct_lbl))
  boolnet_fname<- file.path(out_dir_boolnet, sprintf("ModifiedToyModel_%s.txt",  pct_lbl))
  
  # 3) Save RData & SIF
  save(mod_model, file=rdata_fname)
  writeSIF(mod_model, file=sif_fname, overwrite=TRUE)
  message("Wrote:\n - RData → ", rdata_fname, "\n - SIF   → ", sif_fname, "\n")
  
  # 4) Convert to BoolNet format
  SIFToBoolNet(sifFile     = sif_fname,
               boolnetFile = boolnet_fname,
               CNOlist     = CNOlist,
               model       = orig_model,
               fixInputs   = TRUE,
               preprocess  = TRUE,
               ignoreAnds  = TRUE)
  message("BoolNet file written to: ", boolnet_fname)
  
  # Return paths & perturbed model
  invisible(list(mod_model     = mod_model,
                 sif_fname     = sif_fname,
                 rdata_fname   = rdata_fname,
                 boolnet_fname = boolnet_fname))
}


results <- runPerturbPipeline(
#                 orig_model      = myCustomModel,
#                 CNOlist         = myCNOlist,
                change_pct      = 0.3,
#                 seed            = 123,
#                 out_dir_sif     = "myoutputs/sif",
#                 out_dir_boolnet = "myoutputs/boolnet"
)

print(results)