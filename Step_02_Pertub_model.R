#!/usr/bin/env Rscript
library(CellNOptR)
library(BoolNet)
library(here)
source(here::here("tools", "functions.R"))

#' Perturb a Boolean-logic model by random local edits
#' 
#' @param model  A CellNOptR model list (interMat, notMat, namesSpecies, reacID)
#' @param change_percent  A number in [0,1].  Fraction of nodes to perturb.
#' @return  A modified copy of `model`

perturbModel <- function(model, cnolist, change_percent=0.1, DELETE=FALSE) {
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
    if (DELETE) {
      insertPr <- 1/3
      deletePr <- 2/3
      exchangePr <- 1
    } else {
      insertPr <- 0.5
      deletePr <- 0
      exchangePr <- 1
    }

    ## 1) Operator insertion
    if (x < insertPr) {
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
      
    ## 2) Subtree deletion
    } else if (x < deletePr) {
      message(sprintf("Deleting subtree for node %s", node1))
      # find all reactions involving node1
      cols <- which(grepl(node1, model$reacID, fixed=TRUE))
      message(sprintf("Found %d reactions [%s] for node %s", length(cols), paste(model$reacID[cols], collapse=", "), node1))
      if (length(cols) > 1) {
          col <- sample(cols, 1)  # pick one at random
          message(sprintf("Deleting reaction %s", model$reacID[col]))
          rid <- model$reacID[col]
          model <- remove_reaction(model, col)
      }
      
    ## 3) Operator exchange
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
#' @param dataset      Name of the dataset to use. Options: "toy", "apoptosis", "dream", "TCell".
#'                     Defaults to "toy".
#' @param change_pct   Fraction in [0,1] of nodes to perturb. Default: 0.2.
#' @param seed         Integer RNG seed for reproducibility. Default: 42.
#' @param DELETE       Logical, if TRUE, allows deletion perturbations. Default: FALSE.
#' @return             A list with elements `mod_model`, `sif_fname`, `rdata_fname`, `boolnet_fname`.
#------------------------------------------------------------------------------#
runPerturbPipeline <- function(dataset     = "toy",
                               change_pct   = 0.9,
                               seed         = 42,
                               DELETE       = FALSE,
                               K_FOLD       = 10) {
  # Set the random seed for reproducible results
  # This ensures your cross-validation splits are consistent across runs
  set.seed(seed)
  # 1) Setup dataset mapping (same as your original code)
  # This maps dataset names to their corresponding file names
  dataset_map <- list(
    toy       = c("ToyModel", "ToyModel.sif", "ToyModel.RData", "ToyModel.csv", "ToyModel.bnet"),
    apoptosis = c("Apoptosis", "Apoptosis.sif", "Apoptosis.RData", "Apoptosis.csv", "Apoptosis.bnet"),
    dream     = c("DREAMmodel", "DreamModel.sif", "DreamModel.RData", "DreamModel.csv", "DreamModel.bnet"),
    TCell     = c("T-Cell", "TCell.sif", "TCell.RData", "TCell.csv", "TCell.bnet")
  )
  if (!dataset %in% names(dataset_map)) {
    stop(sprintf("Unknown dataset: %s", dataset))
  }
  
  # Extract file names for the specified dataset
  vals <- dataset_map[[dataset]]
  base_name <- vals[1] 
  sif_name <- vals[2]
  rdata_name <- vals[3]
  midas_name <- vals[4]
  boolnet_name <- vals[5]

  # 2) Create main output directory structure
  # This creates a directory based on the percentage of modification
  new_d  <- paste0(round(change_pct * 100), "_Modified")
  output_file <- file.path("data", base_name, new_d)
  if (!dir.exists(output_file)) dir.create(output_file, recursive=TRUE)
  
  # Define file paths for original and modified files
  sif_fname     <- file.path(output_file, sif_name)
  rdata_fname   <- file.path(output_file, rdata_name)
  midas_fname   <- file.path(output_file, midas_name)
  boolnet_fname <- file.path(output_file, boolnet_name)
  
  # Original file paths
  GD_SIF   <- file.path("data", base_name, sif_name)
  GD_DATA  <- file.path("data", base_name, rdata_name)
  GD_MIDAS <- file.path("data", base_name, midas_name)
  GD_BNET  <- file.path("data", base_name, boolnet_name)

  # Copy the original MIDAS file to the output directory
  file.copy(GD_MIDAS, midas_fname, overwrite=TRUE)

  # 3) Load and perturb the network model
  # Read the original network structure from SIF format
  message("Loading original network model...")
  orig_model <- readSIF(GD_SIF)
  
  # Create CNOlist object from MIDAS data
  # This contains the experimental conditions and measurements
  cnolist <- makeCNOlist(readMIDAS(GD_MIDAS, verbose=TRUE), subfield=FALSE)
  
  # Convert to BoolNet format if needed
  if (!file.exists(GD_BNET)) {
    message("Converting SIF to BoolNet format...")
    # SIFToBoolNet(sifFile     = GD_SIF,
    #              boolnetFile = GD_BNET,
    #              CNOlist     = cnolist,
    #              model       = orig_model,
    #              fixInputs   = FALSE,
    #              preprocess  = TRUE,
    #              ignoreAnds  = TRUE)
    # Convert the model
    result <- writeBnetFromModel(orig_model, GD_BNET) # May have problem

    # Verify the conversion
    verifyBoolNetConversion(orig_model, GD_BNET)
  }

  # Apply perturbations to create modified model
  message("Applying network perturbations...")
  cat(">> DEBUG: number of nodes available = ", length(orig_model$namesSpecies), "\n")
  cat(">> DEBUG: change_pct = ", change_pct, "\n")
  wanted_size <- floor(change_pct * length(orig_model$namesSpecies))
  cat(">> DEBUG: computed size = ", wanted_size, "\n")
  
  # This is where network gets modified according to change_pct
  mod_model <- perturbModel(orig_model, cnolist, change_pct, DELETE=DELETE)
  message(sprintf("Reactions before: %d", length(orig_model$reacID)))
  message(sprintf("Reactions after : %d\n",  length(mod_model$reacID)))
  
  # 4) Save the modified model files
  message("Saving modified model files...")
  save(mod_model, file=rdata_fname)
  writeSIF(mod_model, file=sif_fname, overwrite=TRUE)
  message("Wrote:\n - RData → ", rdata_fname, "\n - SIF   → ", sif_fname, "\n")
  
  # SIFToBoolNet(sifFile     = sif_fname,
              #  boolnetFile = boolnet_fname,
              #  CNOlist     = cnolist,
              #  model       = mod_model,
              #  fixInputs   = FALSE,
              #  preprocess  = TRUE,
              #  ignoreAnds  = TRUE)

  # Convert modified model to BoolNet format
  result <- writeBnetFromModel(mod_model, boolnet_fname)

  # Verify the conversion
  verifyBoolNetConversion(mod_model, boolnet_fname)
  message("BoolNet file written to: ", boolnet_fname)
  
  message("Modified model files saved to: ", output_file)

  # 5) Prepare for k-fold cross-validation
  message("Setting up k-fold cross-validation...")
  
  # Clean up model names (replace hyphens with underscores for compatibility)
  mod_model$reacID <- gsub("-", "_", mod_model$reacID, fixed = TRUE)
  mod_model$namesSpecies <- gsub("-", "_", mod_model$namesSpecies, fixed = TRUE)
  rownames(mod_model$interMat) <- gsub("-", "_", rownames(mod_model$interMat), fixed = TRUE)
  colnames(mod_model$interMat) <- gsub("-", "_", colnames(mod_model$interMat), fixed = TRUE)
  rownames(mod_model$notMat) <- gsub("-", "_", rownames(mod_model$notMat), fixed = TRUE)
  colnames(mod_model$notMat) <- gsub("-", "_", colnames(mod_model$notMat), fixed = TRUE)
  
  # Get the number of experimental samples
  numSamples <- nrow(cnolist@cues)
  message(sprintf("Total samples for cross-validation: %d", numSamples))
  message(sprintf("Creating %d-fold splits with %d repetition(s)", nfold, ntimes))
  
  # 6) Create cross-validation splits
  # This creates random partitions of your data for each run
  set.seed(seed)
  splits <- lapply(1:ntimes, function(run) {
    # Create a random permutation of sample indices
    permut <- sample(1:numSamples, numSamples, replace = FALSE)
    
    # Divide samples into nfold groups
    # Each fold gets approximately equal number of samples
    indices <- lapply(1:nfold, function(i) {
      permut[seq(i, numSamples, nfold)]
    })
    return(indices)
  })
  
  # 7) Create subdirectories for each cross-validation fold
  message("Creating cross-validation subdirectories...")
  for (i in 1:ntimes) {
    for (j in 1:nfold) {
      # Create directory structure: cv_1, cv_2, cv_3, etc.
      cv_dir <- file.path(output_file, paste0("cv_", j))
      if (!dir.exists(cv_dir)) {
        dir.create(cv_dir, recursive = TRUE)
      }
      
      # Save seed information for reproducibility
      seed_file <- file.path(cv_dir, "seed.txt")
      cat(seed, file = seed_file)
    }
  }
  
  # 8) Create train/validation splits for each fold
  message("Creating training and validation datasets for each fold...")
  
  for (i in 1:ntimes) {
    for (j in 1:nfold) {
      cv_dir <- file.path(output_file, paste0("cv_", j))
      
      # Get test indices for this fold (validation set)
      testIndices <- splits[[i]][[j]]
      # Training set includes all other samples
      trainIndices <- (1:numSamples)[-testIndices]
      
      message(sprintf("Fold %d: %d training samples, %d validation samples", 
                     j, length(trainIndices), length(testIndices)))
      
      # Create separate CNOlist objects for training and validation
      CNOlistTrain <- cnolist
      CNOlistVal <- cnolist
      
      # Split the experimental conditions (cues, stimuli, inhibitors)
      CNOlistTrain@cues <- CNOlistTrain@cues[trainIndices, , drop = FALSE]
      CNOlistTrain@stimuli <- CNOlistTrain@stimuli[trainIndices, , drop = FALSE]
      CNOlistTrain@inhibitors <- CNOlistTrain@inhibitors[trainIndices, , drop = FALSE]
      
      CNOlistVal@cues <- CNOlistVal@cues[testIndices, , drop = FALSE]
      CNOlistVal@stimuli <- CNOlistVal@stimuli[testIndices, , drop = FALSE]
      CNOlistVal@inhibitors <- CNOlistVal@inhibitors[testIndices, , drop = FALSE]

      # Split the measurement signals for each time point
      for (k in 1:length(CNOlistTrain@signals)) {
        CNOlistTrain@signals[[k]] <- CNOlistTrain@signals[[k]][trainIndices, , drop = FALSE]
        CNOlistVal@signals[[k]] <- CNOlistVal@signals[[k]][testIndices, , drop = FALSE]
      }
      
      # Save the split datasets
      train_file <- file.path(cv_dir, "CNOlist_train.RData")
      val_file <- file.path(cv_dir, "CNOlist_val.RData")
      model_file <- file.path(cv_dir, "model.RData")
      
      save(CNOlistTrain, file = train_file)
      save(CNOlistVal, file = val_file)

      message(sprintf("Saved fold %d data to: %s", j, cv_dir))
    }
  }
  
  # 9) Save the split information for later reference
  splits_file <- file.path(output_file, "splits.RData")
  save(splits, file = splits_file)
  
  # Create a summary of the cross-validation setup
  cv_summary <- list(
    dataset = dataset,
    change_pct = change_pct,
    nfold = nfold,
    ntimes = ntimes,
    seed = seed,
    numSamples = numSamples,
    output_dir = output_file,
    splits = splits
  )
  
  summary_file <- file.path(output_file, "cv_summary.RData")
  save(cv_summary, file = summary_file)
  
  message("Cross-validation setup complete!")
  message(sprintf("Main directory: %s", output_file))
  message(sprintf("Created %d cross-validation folds", nfold))
  message(sprintf("Each fold contains training and validation datasets"))
  
  return(cv_summary)
}

#--- Load required libraries -------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
})

#--- Define command-line options ---------------------------------------------
option_list <- list(
  make_option(c("-d", "--dataset"), type="character", default="toy",
              help="Path to original model RData or RDS file", metavar="FILE"),
  make_option(c("-p", "--changePCT"), type="double", default=0.9,
              help="Change percentage [default %default]", metavar="DOUBLE"),
  make_option(c("-D", "--delete"), type="logical", default=FALSE,
              help="Enable deletion perturbation [default %default]", metavar="LOGICAL"),
  make_option(c("-s", "--seed"), type="integer", default=44,
              help="Random seed [optional]", metavar="INT"),
  make_option(c("-k", "--k_fold"), type="integer", default=10,
              help="Number of folds for cross-validation [default %default]", metavar="INT")
)

parser <- OptionParser(option_list=option_list,
                       description = "Run perturbation pipeline for a Boolean model")
opt <- parse_args(parser)

#--- (Optional) set seed -----------------------------------------------------
if (!is.null(opt$seed)) {
  set.seed(opt$seed)
  message("Random seed set to: ", opt$seed)
}

#--- Run the pipeline --------------------------------------------------------
results <- runPerturbPipeline(
    dataset      = opt$dataset,
    change_pct   = opt$changePCT,
    DELETE       = opt[["delete"]],
    seed         = opt$seed
  )
