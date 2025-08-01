#!/usr/bin/env Rscript
library(CellNOptR)
library(BoolNet)
library(here)
source(here::here("tools", "functions.R"))

#===============================================================================
# OPTIMIZED NETWORK PERTURBATION SYSTEM FOR PKN MODIFICATION
#===============================================================================

#' Core function to apply intelligent perturbations to a biological network model
#' 
#' This function implements sophisticated perturbation strategies that test different
#' types of biological uncertainty, going beyond simple edge removal to include
#' sign flipping, direction changes, pathway rewiring, and relationship splitting/merging.
#'
#' @param orig_model Original CNO model object with interMat, notMat, reacID, namesSpecies
#' @param cnolist CNOlist object containing experimental conditions  
#' @param change_pct Fraction of relationships to modify (0-1)
#' @param seed Random seed for reproducibility
#' @return Modified model object with perturbed network structure
perturbModel <- function(orig_model, cnolist, change_pct = 0.2, seed = 42) {
  set.seed(seed)  # Set the random seed for reproducibility
  message(sprintf("Setting random seed to %d for reproducibility", seed))
  # Set the random seed for reproducibility  # Create working copy of the model
  mod_model <- orig_model
  
  # Extract key network components for perturbation logic
  network_info <- extractNetworkInfo(mod_model, cnolist)
  n_relations  <- length(mod_model$reacID)
  
  # Extract relationships to modify based on the specified change percentage
  candidates   <- extractModifiableRelations(mod_model, network_info)
  n_candidates <- length(candidates)
  n_to_modify  <- max(1, floor(change_pct * n_candidates))
  
  message(sprintf(
    "Modifying %d out of %d candidate relationships (%.1f%% of %d total)",
    n_to_modify, n_candidates, change_pct * 100, n_relations
  ))

  # Apply perturbations to each selected relationship
  for (i in seq_len(n_to_modify)) {
    candidates <- extractModifiableRelations(mod_model, network_info)
    rel <- sample(candidates, size = 1)
    rel_idx <- match(rel, mod_model$reacID)
    message(sprintf(
      "[%d/%d] Perturbing reaction '%s'/'%s'",
      i, n_to_modify, rel_idx, rel
    ))
    mod_model <- applyPerturbation(mod_model, rel_idx, network_info)
  }
  
  # Ensure model integrity after all modifications
  mod_model <- validateModel(mod_model)
  
  return(mod_model)
}

# ' Extract relationships that can be modified based on biological context
extractModifiableRelations <- function(model, network_info, change_pct) {
  # For each reaction ID, find what perturbations are possible
  rel_list <- lapply(seq_along(model$reacID), function(idx) {
    rxn_id <- model$reacID[idx]
    info   <- parseRelationship(model, idx)
    ap     <- getAvailablePerturbations(model, info, network_info)
    
    message(sprintf("Available perturbations for %s: %s", 
                    rxn_id, if (length(ap)) paste(ap, collapse = ", ") else "<none>"))
    if (length(ap) > 0) {
      rxn_id
    } else {
      NULL
    }
  })
  
  # Drop the NULLs, return a character vector of IDs
  Filter(Negate(is.null), rel_list)
}

#' Extract key network information needed for intelligent perturbation decisions
#' 
#' This helper function analyzes the network structure to identify stimuli, inhibitors,
#' readouts, and leaf nodes, which are used to make context-aware perturbation choices.
extractNetworkInfo <- function(model, cnolist) {
  # Get node types from experimental setup
  stimuli <- cnolist$namesStimuli
  inhibitors <- cnolist$namesInhibitors
  readouts <- cnolist$namesSignals
  
  # Clean up names to match model species (remove special characters)
  stimuli <- intersect(stimuli, model$namesSpecies)
  inhibitors <- intersect(inhibitors, model$namesSpecies)
  readouts <- intersect(readouts, model$namesSpecies)
  
  # Identify leaf nodes (nodes that are only targets, never sources)
  all_sources <- getSourceNodes(model)
  all_targets <- getTargetNodes(model)
  leaf_nodes <- setdiff(all_targets, all_sources)
  
  # Create comprehensive network info structure
  list(
    stimuli = stimuli,
    inhibitors = inhibitors, 
    readouts = readouts,
    leaf_nodes = leaf_nodes,
    all_nodes = model$namesSpecies,
    all_sources = all_sources,
    all_targets = all_targets
  )
}

#' Apply a single perturbation to a relationship based on biological context
#' 
#' This function chooses the appropriate perturbation type based on the biological
#' role of the source node (stimulus, inhibitor, readout, or regular node).
applyPerturbation <- function(model, rel_idx, network_info) {
  # Parse the current relationship
  rel_info <- parseRelationship(model, rel_idx)
  
  # Get available perturbation types based on biological context
  available_perturbations <- getAvailablePerturbations(model, rel_info, network_info)
  
  # Randomly select a perturbation type
  perturbation_type <- sample(available_perturbations, 1)
  
  # Apply the selected perturbation
  model <- executePerturbation(model, rel_idx, rel_info, perturbation_type, network_info)
  
  return(model)
}

#' Determine which perturbation types are allowed based on biological context
#' 
#' Different types of nodes (stimuli, inhibitors, readouts) have different constraints
#' on how they can be perturbed while maintaining biological plausibility.
getAvailablePerturbations <- function(model, rel_info, network_info) {
  base_perturbations <- c("flip_sign", "flip_direction", "flip_sign_and_direction")
  
  # Add complex perturbations if the network structure supports them
  if (canChangeTarget(model, rel_info, network_info)) {
    base_perturbations <- c(base_perturbations, "change_target")
  }
  if (canChangeSource(model, rel_info, network_info)) {
    base_perturbations <- c(base_perturbations, "change_source")
  }
  if (canRewirePathway(model, rel_info, network_info)) {
    base_perturbations <- c(base_perturbations, "rewire_pathway")
  }
  if (canSplitRelationship(rel_info, network_info)) {
    base_perturbations <- c(base_perturbations, "split_relationship")
  }  
  # if (canMergeRelationship(model, rel_info, network_info)) {
  #   base_perturbations <- c(base_perturbations, "merge_relationships")
  # }
  
  # Context-specific filtering based on node types
  if (rel_info$source %in% network_info$stimuli) {
    # Stimuli nodes: can only change targets, not sources
    return(intersect(base_perturbations, c("change_target")))
    
  } else if (rel_info$source %in% network_info$inhibitors) {
    # Inhibitor nodes: limited perturbation options    
    forbidden <- c() 
    return(setdiff(base_perturbations, forbidden))

  } else if (rel_info$target %in% network_info$leaf_nodes) {
    # Leaf nodes: cannot change source (would create dangling edges)
    forbidden <- c() 
    return(setdiff(base_perturbations, forbidden))
    
  } else {
    # Regular internal nodes: all perturbations allowed
    return(base_perturbations)
  }
}

#' Execute the actual perturbation operation on the model
#' 
#' This is the main dispatcher that calls the appropriate modification function
#' based on the perturbation type selected.
executePerturbation <- function(model, rel_idx, rel_info, perturbation_type, network_info) {
  message(sprintf("Executing perturbation: %s on relationship %s", 
                  perturbation_type, rel_idx))
  switch(perturbation_type,
    "change_target" = changeTarget(model, rel_idx, rel_info, network_info),
    "change_source" = changeSource(model, rel_idx, rel_info, network_info),
    "flip_sign" = flipSign(model, rel_idx, rel_info),
    "flip_direction" = flipDirection(model, rel_idx, rel_info),
    "flip_sign_and_direction" = flipSignAndDirection(model, rel_idx, rel_info),
    "rewire_pathway" = rewirePathway(model, rel_idx, rel_info, network_info),
    "split_relationship" = splitRelationship(model, rel_idx, rel_info, network_info),
    "merge_relationships" = mergeRelationships(model, rel_idx, rel_info, network_info),
    model  # Return unchanged if perturbation type not recognized
  )
}

#===============================================================================
# SPECIFIC PERTURBATION IMPLEMENTATIONS
#===============================================================================

#' Change the target of a relationship while keeping the same source
#' 
#' This tests uncertainty about which genes are regulated by a given regulator.
changeTarget <- function(model, rel_idx, rel_info, network_info) {
  # Find connected targets to avoid creating duplicate relationships
  connected_targets <- getConnectedTargets(model, rel_info$source)
  
  # Select new target from available nodes (excluding stimuli and current connections)
  possible_targets <- setdiff(network_info$all_nodes, 
                             c(network_info$stimuli, connected_targets, rel_info$source))
    
  new_target <- sample(possible_targets, 1)
  
  # Remove old relationship and add new one
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, rel_info$source, new_target, rel_info$is_inhibition)
  
  return(model)
}

#' Change the source of a relationship while keeping the same target
#' 
#' This tests uncertainty about which genes regulate a given target.
changeSource <- function(model, rel_idx, rel_info, network_info) {
  # Find connected sources to avoid creating duplicate relationships  
  connected_sources <- getConnectedSources(model, rel_info$target)
  
  # Select new source from available nodes (excluding stimuli and current connections)
  possible_sources <- setdiff(network_info$all_nodes,
                             c(network_info$stimuli, connected_sources, rel_info$target))
  
  new_source <- sample(possible_sources, 1)
  
  # Remove old relationship and add new one
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, new_source, rel_info$target, rel_info$is_inhibition)
  
  return(model)
}

#' Flip the sign of a relationship (activation ↔ inhibition)
#' 
#' This tests uncertainty about whether interactions are activating or inhibiting.
flipSign <- function(model, rel_idx, rel_info) {
  # Simply toggle the inhibition status
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, rel_info$source, rel_info$target, !rel_info$is_inhibition)
  
  return(model)
}

#' Flip the direction of a relationship (A→B becomes B→A)
#' 
#' This tests uncertainty about causality and feedback relationships.
flipDirection <- function(model, rel_idx, rel_info) {
  # Swap source and target
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, rel_info$target, rel_info$source, rel_info$is_inhibition)
  
  return(model)
}

#' Flip both sign and direction of a relationship
#' 
#' This combines sign and direction uncertainty for maximum perturbation.
flipSignAndDirection <- function(model, rel_idx, rel_info) {
  # Swap source/target and toggle inhibition
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, rel_info$target, rel_info$source, !rel_info$is_inhibition)
  
  return(model)
}

#' Rewire pathway connections to test pathway order sensitivity
#' 
#' This finds a pathway like A→B→C and changes it to A→C→B to test
#' how sensitive methods are to the specific order of regulatory events.
rewirePathway <- function(model, rel_idx, rel_info, network_info) {
  # Find a downstream relationship from current target
  downstream_rels <- findDownstreamRelationships(model, rel_info$target, k = 1)
  
  # Select one downstream relationship to rewire
  model <- removeReaction(model, rel_idx)

  downstream_rel <- sample(downstream_rels, 1)
  downstream_idx <- match(downstream_rel, model$reacID)
  message(sprintf("Rewire downstream relationship: %s", downstream_idx))
  downstream_info <- parseRelationship(model, downstream_idx)
  message(sprintf("Rewire %s with %s", rel_info$reaction_id, downstream_info$reaction_id))

  # Create the rewired connections: A→C and C→B instead of A→B→C
  model <- removeReaction(model, downstream_idx)
  
  # Add new connections with original signs
  model <- addReaction(model, rel_info$source, downstream_info$target, rel_info$is_inhibition)
  model <- addReaction(model, downstream_info$target, rel_info$target, downstream_info$is_inhibition)
  
  return(model)
}

#' Split a direct relationship into an indirect one through an intermediate node
#' 
#' This replaces A→C with A→B→C to test how methods handle direct versus indirect regulation.
splitRelationship <- function(model, rel_idx, rel_info, network_info) {
  # Find a suitable intermediate node
  potential_intermediates <- setdiff(network_info$all_nodes, 
                                   c(rel_info$source, rel_info$target, 
                                     network_info$stimuli, network_info$leaf_nodes))
  
  intermediate <- sample(potential_intermediates, 1)
  
  # Remove direct relationship and add indirect path
  model <- removeReaction(model, rel_idx)
  model <- addReaction(model, rel_info$source, intermediate, rel_info$is_inhibition)
  model <- addReaction(model, intermediate, rel_info$target, FALSE)  # Second step as activation
  
  return(model)
}

#' Merge two sequential relationships into one direct relationship
#' 
#' This finds A→B→C and merges it into A→C to test how methods handle
#' merging multiple regulatory relationships into a single one.
mergeRelationships <- function(model, rel_idx, rel_info, network_info) {
  # Look for relationships that could be merged with the current one
  if (rel_info$target %in% network_info$leaf_nodes) {
    # Current target is a leaf, look for upstream relationships
    upstream_rels <- findUpstreamRelationships(model, rel_info$source, k = 1)    

    # Remove both relationships and create merged one
    model <- removeReaction(model, rel_idx)

    merge_rel <- sample(upstream_rels, 1)
    merge_idx <- match(merge_rel, model$reacID)
    message(sprintf("Merging upstream relationship: %s", merge_idx))
    merge_info <- parseRelationship(model, merge_idx)
    message(sprintf("Merging %s with %s", rel_info$reaction_id, merge_info$reaction_id))
    
    model <- removeReaction(model, merge_idx)
    model <- addReaction(model, merge_info$source, rel_info$target, 
                        merge_info$is_inhibition)  # Keep upstream sign
    
  } else {
    # Look for downstream relationships to merge
    downstream_rels <- findDownstreamRelationships(model, rel_info$target, k = 1)
    
    # Remove both relationships and create merged one
    model <- removeReaction(model, rel_idx)

    merge_rel <- sample(downstream_rels, 1)
    merge_idx <- match(merge_rel, model$reacID)
    message(sprintf("Merging upstream relationship: %s", merge_idx))
    merge_info <- parseRelationship(model, merge_idx)
    message(sprintf("Merging %s with %s", rel_info$reaction_id, merge_info$reaction_id))
    
    model <- removeReaction(model, merge_idx)
    model <- addReaction(model, rel_info$source, merge_info$target, 
                        rel_info$is_inhibition)  # Keep original sign
  }
  
  return(model)
}

#===============================================================================
# HELPER FUNCTIONS FOR MODEL MANIPULATION
#===============================================================================

#' Parse a relationship from the model into a structured format
parseRelationship <- function(model, rel_idx) {
  reaction_id <- model$reacID[rel_idx]
  
  # Handle negation prefix (!)
  is_inhibition <- startsWith(reaction_id, "!")
  clean_id <- sub("^!", "", reaction_id)
  
  # Split on = to get source and target
  parts <- strsplit(clean_id, "=", fixed = TRUE)[[1]]
  
  list(
    source = parts[1],
    target = parts[2], 
    is_inhibition = is_inhibition,
    reaction_id = reaction_id
  )
}

#' Add a new reaction to the model with proper matrix updates
addReaction <- function(model, source, target, is_inhibition = FALSE) {
  message(sprintf("Adding reaction: %s → %s (Inhibition: %s)", 
                  source, target, ifelse(is_inhibition, "Yes", "No")))
  # Create new reaction ID
  new_id <- if (is_inhibition) paste0("!", source, "=", target) else paste0(source, "=", target)
  
  # Skip if reaction already exists
  if (new_id %in% model$reacID) return(model)
  
  # Add to reaction ID list
  model$reacID <- c(model$reacID, new_id)
  
  # Extend matrices with new column
  n_species <- length(model$namesSpecies)
  model$interMat <- cbind(model$interMat, rep(0, n_species))
  model$notMat <- cbind(model$notMat, rep(0, n_species))
  
  # Set matrix values for new reaction
  col_idx <- ncol(model$interMat)
  source_idx <- match(source, model$namesSpecies)
  target_idx <- match(target, model$namesSpecies)
  
  model$interMat[source_idx, col_idx] <- -1  # Source
  model$interMat[target_idx, col_idx] <- +1  # Target
  
  if (is_inhibition) {
    model$notMat[source_idx, col_idx] <- 1  # Mark as inhibition
  }
  
  # Update column names
  colnames(model$interMat) <- model$reacID
  colnames(model$notMat) <- model$reacID
  
  return(model)
}

#' Remove a reaction from the model by column index
removeReaction <- function(model, col_idx) {
  message(sprintf("Removing reaction %s", model$reacID[col_idx]))
  model$reacID <- model$reacID[-col_idx]
  model$interMat <- model$interMat[, -col_idx]
  model$notMat <- model$notMat[, -col_idx]
  
  # Update column names
  colnames(model$interMat) <- model$reacID
  colnames(model$notMat) <- model$reacID
  
  return(model)
}

#' Get all source nodes in the network
getSourceNodes <- function(model) {
  source_indices <- which(model$interMat == -1, arr.ind = TRUE)[, 1]
  unique(model$namesSpecies[source_indices])
}

#' Get all target nodes in the network  
getTargetNodes <- function(model) {
  target_indices <- which(model$interMat == +1, arr.ind = TRUE)[, 1]
  unique(model$namesSpecies[target_indices])
}

#' Find all targets connected to a given source
getConnectedTargets <- function(model, source) {
  source_idx <- match(source, model$namesSpecies)
  reaction_cols <- which(model$interMat[source_idx, ] == -1)
  target_indices <- apply(model$interMat[, reaction_cols, drop = FALSE], 2, 
                         function(col) which(col == +1))
  unique(model$namesSpecies[unlist(target_indices)])
}

#' Find all sources connected to a given target
getConnectedSources <- function(model, target) {
  target_idx <- match(target, model$namesSpecies)
  reaction_cols <- which(model$interMat[target_idx, ] == +1)
  source_indices <- apply(model$interMat[, reaction_cols, drop = FALSE], 2,
                         function(col) which(col == -1))
  unique(model$namesSpecies[unlist(source_indices)])
}

#' Check if a relationship can change target (has suitable sources)
canChangeTarget <- function(model, rel_info, network_info) {
  connected_sources <- getConnectedSources(model, rel_info$target)
  
  # Select new source from available nodes (excluding stimuli and current connections)
  possible_sources <- setdiff(network_info$all_nodes,
                             c(network_info$stimuli, connected_sources, rel_info$target))

  return(length(possible_sources) > 0 && length(connected_sources) > 1)
}

#' Check if a relationship can change source (has suitable targets)
canChangeSource <- function(model, rel_info, network_info) {
  connected_targets <- getConnectedTargets(model, rel_info$source)
  
  # Select new target from available nodes (excluding stimuli and current connections)
  possible_targets <- setdiff(network_info$all_nodes, 
                             c(network_info$stimuli, connected_targets, rel_info$source))

  return(length(possible_targets) > 0 && length(connected_targets) > 1)
}

#' Check if a relationship can be split (has suitable intermediate nodes)
canSplitRelationship <- function(rel_info, network_info) {
  potential_intermediates <- setdiff(network_info$all_nodes,
                                   c(rel_info$source, rel_info$target,
                                     network_info$stimuli, network_info$leaf_nodes))
  length(potential_intermediates) > 0
}

#' Check if a relationship can be rewired (has downstream relationships)
canRewirePathway <- function(model, rel_info, network_info) {
  # Check if there are downstream relationships to rewire
  downstream_rels <- findDownstreamRelationships(model, rel_info$target, k = 1)

  has_downstream <- length(downstream_rels) > 0

  return(has_downstream)
}
#' Check if a relationship can be merged (has connected relationships)
canMergeRelationship <- function(model, rel_info, network_info) {
  source_node <- rel_info$source
  target_node <- rel_info$target
  
  if (target_node %in% network_info$leaf_nodes) {
    # Current target is a leaf, look for upstream relationships
    upstream_rels <- findUpstreamRelationships(model, source_node, k = 1)    
    
    has_upstream <- length(upstream_rels) > 0

    return(has_upstream)
  } else {
  
    # Run k-hop search in both directions
    downstream_rels <- findDownstreamRelationships(model, target_node, k = 1)

    # Only allow merging if both k-hop neighbors exist
    has_downstream <- length(downstream_rels) > 0

    return(has_downstream)
  }
}

#' Find downstream relationships from a given node
findDownstreamRelationships <- function(model, node, k=1) {
  if (k <= 0) return(character(0))
  
  # find all reactions where 'node' is the source
  rel_idxs <- which(sapply(seq_along(model$reacID), function(i) {
    parseRelationship(model, i)$source == node
  }))
  
  # for each such reaction, record it and recurse on its target
  out <- unlist(lapply(rel_idxs, function(i) {
    rel <- parseRelationship(model, i)
    c(rel$reaction_id,
      findDownstreamRelationships(model, rel$target, k - 1))
  }), use.names = FALSE)
  
  unique(out)
}

#' Find upstream relationships to a given node
findUpstreamRelationships <- function(model, node, k=1) {
  if (k <= 0) return(character(0))
  
  # find all reactions where 'node' is the target
  rel_idxs <- which(sapply(seq_along(model$reacID), function(i) {
    parseRelationship(model, i)$target == node
  }))
  
  # for each such reaction, record it and recurse on its source
  out <- unlist(lapply(rel_idxs, function(i) {
    rel <- parseRelationship(model, i)
    c(rel$reaction_id,
      findUpstreamRelationships(model, rel$source, k - 1))
  }), use.names = FALSE)
  
  unique(out)
}

#' Validate model integrity after modifications
validateModel <- function(model) {
  # Remove any empty reactions or invalid references
  valid_reactions <- !is.na(model$reacID) & model$reacID != ""
  
  model$reacID <- model$reacID[valid_reactions]
  model$interMat <- model$interMat[, valid_reactions, drop = FALSE]
  model$notMat <- model$notMat[, valid_reactions, drop = FALSE]
  
  # Update column names
  colnames(model$interMat) <- model$reacID
  colnames(model$notMat) <- model$reacID
  
  return(model)
}


#------------------------------------------------------------------------------#
#   Function: runPerturbPipeline                                                #
#------------------------------------------------------------------------------#
#' @param dataset      Name of the dataset to use. Options: "toy", "apoptosis", "dream", "TCell".
#'                     Defaults to "toy".
#' @param change_pct   Fraction in [0,1] of nodes to perturb. Default: 0.2.
#' @param seed         Integer RNG seed for reproducibility. Default: 42.
#' @param nfold        Number of folds for cross-validation. Default: 10.
#' @return             A list with elements `mod_model`, `sif_fname`, `rdata_fname`, `boolnet_fname`.
#------------------------------------------------------------------------------#
runPerturbPipeline <- function(dataset      = "toy",
                               change_pct   = 0.9,
                               seed         = 42,
                               nfold        = 10) {
  # Set the random seed for reproducible results
  # This ensures your cross-validation splits are consistent across runs

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
  cat(">> DEBUG: using seed = ", seed, "\n")
  wanted_size <- floor(change_pct * length(orig_model$namesSpecies))
  cat(">> DEBUG: computed size = ", wanted_size, "\n")
  
  # This is where network gets modified according to change_pct
  set.seed(seed)
  mod_model <- perturbModel(orig_model, cnolist, change_pct, seed)
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
  numSamples <- length(cnolist$namesCues)
  nfold <- min(nfold, numSamples)  # Ensure nfold does not exceed number of samples
  message(sprintf("Total samples for cross-validation: %d", numSamples))
  message(sprintf("Creating %d-fold splits", nfold))

  # 6) Create cross-validation splits
  # This creates random partitions of your data for each run
  cv_seed <- seed + 1000  # Use different seed for CV splits
  set.seed(cv_seed)
  message(sprintf("Using CV seed: %d", cv_seed))
  # Create a random permutation of sample indices
  permut <- sample(seq_len(numSamples))
  
  # Divide samples into nfold groups
  # assign each position to one of nfold buckets as evenly as possible
  fold_id <- cut(seq_along(permut),
                breaks = nfold,
                labels = FALSE)

  # collect your splits
  splits <- lapply(seq_len(nfold), function(k) {
    permut[fold_id == k]
  })
  
  # 7) Create subdirectories for each cross-validation fold
  message("Creating cross-validation subdirectories...")
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
  
  # 8) Create train/validation splits for each fold
  message("Creating training and validation datasets for each fold...")
  
  for (j in 1:nfold) {
    cv_dir <- file.path(output_file, paste0("cv_", j))
    
    # Get test indices for this fold (validation set)
    testIndices <- splits[[j]]
    # Training set includes all other samples
    trainIndices <- setdiff(seq_len(numSamples), testIndices)

    message(sprintf("Fold %d: %d training samples, %d validation samples",
                    j, length(trainIndices), length(testIndices)))
    
    # Create separate CNOlist objects for training and validation
    CNOlistTrain <- cnolist
    CNOlistVal   <- cnolist

    # Split the experimental conditions (cues, stimuli, inhibitors)
    CNOlistTrain$valueCues       <- cnolist$valueCues[    trainIndices, , drop = FALSE]
    CNOlistVal$valueCues         <- cnolist$valueCues[    testIndices,  , drop = FALSE]

    CNOlistTrain$valueStimuli    <- cnolist$valueStimuli[ trainIndices, , drop = FALSE]
    CNOlistVal$valueStimuli      <- cnolist$valueStimuli[ testIndices,  , drop = FALSE]

    CNOlistTrain$valueInhibitors <- cnolist$valueInhibitors[trainIndices, , drop = FALSE]
    CNOlistVal$valueInhibitors   <- cnolist$valueInhibitors[testIndices,  , drop = FALSE]


    # Split the measurement signals for each time point
    for (k in seq_along(cnolist$valueSignals)) {
      CNOlistTrain$valueSignals[[k]]    <- cnolist$valueSignals[[k]][   trainIndices, , drop = FALSE]
      CNOlistVal$valueSignals[[k]]      <- cnolist$valueSignals[[k]][   testIndices,  , drop = FALSE]
      CNOlistTrain$valueVariances[[k]]  <- cnolist$valueVariances[[k]][ trainIndices, , drop = FALSE]
      CNOlistVal$valueVariances[[k]]    <- cnolist$valueVariances[[k]][ testIndices,  , drop = FALSE]
    }
    
    # Save the split datasets
    train_file <- file.path(cv_dir, "CNOlist_train.RData")
    val_file <- file.path(cv_dir, "CNOlist_val.RData")
    
    save(CNOlistTrain, file = train_file)
    save(CNOlistVal, file = val_file)

    message(sprintf("Saved fold %d data to: %s", j, cv_dir))
  }
  
  # 9) Save the split information for later reference
  splits_file <- file.path(output_file, "splits.RData")
  save(splits, file = splits_file)
  
  # Create a summary of the cross-validation setup
  cv_summary <- list(
    dataset = dataset,
    change_pct = change_pct,
    nfold = nfold,
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
  make_option(c("-s", "--seed"), type="integer", default=44,
              help="Random seed [optional]", metavar="INT"),
  make_option(c("-k", "--k_fold"), type="integer", default=10,
              help="Number of folds for cross-validation [default %default]", metavar="INT")
)

parser <- OptionParser(option_list=option_list,
                       description = "Run perturbation pipeline for a Boolean model")
opt <- parse_args(parser)

message("Running perturbation pipeline with the following parameters:")
message(sprintf("Dataset: %s", opt$dataset))
message(sprintf("Change percentage: %.2f", opt$changePCT))
message(sprintf("Random seed: %d", opt$seed)) 
message(sprintf("Number of folds: %d", opt$k_fold))

#--- Run the pipeline --------------------------------------------------------
results <- runPerturbPipeline(
    dataset      = opt$dataset,
    change_pct   = opt$changePCT,
    seed         = opt$seed,
    nfold        = opt$k_fold
  )
