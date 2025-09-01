library(BoolNet)
library(CellNOptR)
library(here)
library(optparse)
library(igraph)

source(here::here("tools", "fundatactions.R"))

#' Identify node types based on network topology and SBML annotations
#' @param network BoolNet network object or igraph object
#' @return List containing stimuli, inhibitors, readouts, and external nodes
identify_node_types <- function(network) {
  
  # Handle different network types
  if (class(network)[1] == "BooleanNetwork") {
    # Convert BoolNet to igraph for analysis
    adj_matrix <- getInteractionMatrix(network)
    graph <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
  } else if (class(network)[1] == "igraph") {
    graph <- network
  } else {
    stop("Network must be BooleanNetwork or igraph object")
  }
  
  # Get network properties
  node_names <- V(graph)$name
  if (is.null(node_names)) {
    node_names <- as.character(1:vcount(graph))
  }
  
  in_degrees <- degree(graph, mode = "in")
  out_degrees <- degree(graph, mode = "out")
  
  # External nodes (no incoming edges) - these are typically stimuli
  external_nodes <- names(which(in_degrees == 0))
  
  # Stimuli: external nodes + nodes with very low in-degree but high out-degree
  stimuli <- unique(c(
    external_nodes,
    names(which(in_degrees <= 1 & out_degrees >= 2))
  ))
  
  # Readouts: high in-degree, low out-degree (sinks or near-sinks)
  readouts <- names(which(in_degrees >= 2 & out_degrees <= 1))
  
  # Inhibitors: identify from edge signs if available, or use heuristics
  inhibitors <- c()
  if ("weight" %in% edge_attr_names(graph)) {
    # Look for nodes that primarily have negative outgoing edges
    for (node in node_names) {
      out_edges <- E(graph)[from(node)]
      if (length(out_edges) > 0) {
        edge_weights <- edge_attr(graph, "weight", out_edges)
        if (mean(edge_weights < 0) > 0.6) {  # More than 60% negative edges
          inhibitors <- c(inhibitors, node)
        }
      }
    }
  }
  
  # Remove overlap between categories
  readouts <- setdiff(readouts, stimuli)
  inhibitors <- setdiff(inhibitors, c(stimuli, readouts))
  
  # Ensure we have at least some nodes in each category
  if (length(stimuli) == 0) {
    stimuli <- external_nodes
  }
  if (length(readouts) == 0 && length(node_names) > 3) {
    # Pick nodes with highest in-degree as fallback
    top_in_degree <- head(names(sort(in_degrees, decreasing = TRUE)), 
                         max(1, floor(length(node_names) * 0.2)))
    readouts <- setdiff(top_in_degree, stimuli)
  }
  
  return(list(
    stimuli = stimuli,
    inhibitors = inhibitors, 
    readouts = readouts,
    external = external_nodes,
    all_nodes = node_names
  ))
}

#' Design comprehensive experimental conditions
#' @param stimuli Vector of stimulus node names
#' @param inhibitors Vector of inhibitor node names  
#' @param readouts Vector of readout node names
#' @param include_combinations Logical, whether to include combination experiments
#' @return List of experimental conditions
design_experiments <- function(stimuli, inhibitors, readouts, include_combinations = TRUE) {
  experiments <- list()
  all_perturbations <- c(stimuli, inhibitors)
  
  # Control condition (no perturbations)
  control_condition <- rep(0, length(all_perturbations))
  names(control_condition) <- all_perturbations
  experiments[["control"]] <- list(
    perturbations = control_condition,
    readouts = readouts,
    type = "control"
  )
  
  # Single stimulus experiments
  for (stim in stimuli) {
    exp_name <- paste0("stim_", stim)
    condition <- rep(0, length(all_perturbations))
    names(condition) <- all_perturbations
    condition[stim] <- 1
    experiments[[exp_name]] <- list(
      perturbations = condition,
      readouts = readouts,
      type = "stimulus"
    )
  }
  
  # Single inhibitor experiments  
  for (inhib in inhibitors) {
    exp_name <- paste0("inhib_", inhib)
    condition <- rep(0, length(all_perturbations))
    names(condition) <- all_perturbations
    condition[inhib] <- 1  # Inhibitor active (blocking signal)
    experiments[[exp_name]] <- list(
      perturbations = condition,
      readouts = readouts,
      type = "inhibitor"
    )
  }
  
  # Combination experiments (if requested and feasible)
  if (include_combinations && length(stimuli) > 1) {
    # All stimuli combinations (limited to avoid explosion)
    if (length(stimuli) <= 4) {
      for (i in 2:min(length(stimuli), 3)) {
        stim_combos <- combn(stimuli, i, simplify = FALSE)
        for (combo in head(stim_combos, 5)) {  # Limit to 5 combinations
          exp_name <- paste0("combo_", paste(combo, collapse = "_"))
          condition <- rep(0, length(all_perturbations))
          names(condition) <- all_perturbations
          condition[combo] <- 1
          experiments[[exp_name]] <- list(
            perturbations = condition,
            readouts = readouts,
            type = "combination"
          )
        }
      }
    }
  }
  
  return(experiments)
}

#' Perform Markov simulation with proper perturbation handling
#' @param network BoolNet network object
#' @param experiments List of experimental conditions
#' @param n_simulations Number of simulation iterations
#' @param method Simulation method ("random" or "exhaustive")
#' @return Data frame with MIDAS-formatted results
generate_midas_data <- function(network, experiments, n_simulations = 1000, method = "random") {
  
  midas_data <- list()
  
  for (exp_name in names(experiments)) {
    exp <- experiments[[exp_name]]
    
    cat(sprintf("Processing experiment: %s\n", exp_name))
    
    # Create perturbed network by fixing perturbation nodes
    perturbed_network <- network
    
    # Apply perturbations
    perturbation_nodes <- names(exp$perturbations)[exp$perturbations != 0]
    perturbation_values <- exp$perturbations[exp$perturbations != 0]
    
    if (length(perturbation_nodes) > 0) {
      for (i in seq_along(perturbation_nodes)) {
        node <- perturbation_nodes[i]
        value <- perturbation_values[i]
        
        # Check if node exists in network
        if (node %in% names(perturbed_network$genes)) {
          perturbed_network <- fixGenes(perturbed_network, node, value)
        } else {
          warning(sprintf("Node %s not found in network", node))
        }
      }
    }
    
    # Perform simulation based on method
    tryCatch({
      if (method == "exhaustive" && length(perturbed_network$genes) <= 15) {
        # Use attractor analysis for smaller networks
        attractors <- getAttractors(perturbed_network, method = "exhaustive")
        readout_values <- calculate_attractor_readouts(attractors, exp$readouts)
      } else {
        # Use Markov simulation for larger networks
        sim_result <- markovSimulation(
          network = perturbed_network,
          method = "random",
          startStates = min(100, 2^min(length(perturbed_network$genes), 10)),
          numIterations = n_simulations,
          returnTable = TRUE
        )
        
        readout_values <- calculate_markov_readouts(sim_result, exp$readouts)
      }
      
      # Create MIDAS row
      midas_row <- c(
        list(ID = exp_name),
        as.list(exp$perturbations),
        as.list(readout_values)
      )
      
      midas_data[[exp_name]] <- midas_row
      
    }, error = function(e) {
      warning(sprintf("Error in experiment %s: %s", exp_name, e$message))
      # Create row with NA values
      readout_values <- rep(NA, length(exp$readouts))
      names(readout_values) <- exp$readouts
      
      midas_row <- c(
        list(ID = exp_name),
        as.list(exp$perturbations),
        as.list(readout_values)
      )
      
      midas_data[[exp_name]] <- midas_row
    })
  }
  
  # Convert to data frame
  midas_df <- do.call(rbind, lapply(midas_data, function(x) data.frame(x, stringsAsFactors = FALSE)))
  
  return(midas_df)
}

#' Calculate readout values from Markov simulation results
#' @param sim_result Markov simulation result
#' @param readouts Vector of readout node names
#' @return Named vector of readout values
calculate_markov_readouts <- function(sim_result, readouts) {
  steady_states <- sim_result$table
  readout_values <- c()
  
  for (readout in readouts) {
    if (readout %in% colnames(steady_states)) {
      # Weighted average based on state probabilities
      expected_val <- sum(steady_states[, readout] * steady_states$Probability)
      readout_values <- c(readout_values, expected_val)
    } else {
      warning(sprintf("Readout %s not found in simulation results", readout))
      readout_values <- c(readout_values, NA)
    }
  }
  
  names(readout_values) <- readouts
  return(readout_values)
}

#' Calculate readout values from attractor analysis
#' @param attractors Attractor analysis result
#' @param readouts Vector of readout node names
#' @return Named vector of readout values
calculate_attractor_readouts <- function(attractors, readouts) {
  readout_values <- c()
  
  for (readout in readouts) {
    if (readout %in% rownames(attractors$attractors[[1]])) {
      # Average across all attractors weighted by basin size
      total_states <- sum(attractors$basinSize)
      weighted_sum <- 0
      
      for (i in seq_along(attractors$attractors)) {
        attractor <- attractors$attractors[[i]]
        basin_weight <- attractors$basinSize[i] / total_states
        
        if (ncol(attractor) == 1) {
          # Fixed point attractor
          node_val <- attractor[readout, 1]
        } else {
          # Cyclic attractor - take average
          node_val <- mean(attractor[readout, ])
        }
        
        weighted_sum <- weighted_sum + (node_val * basin_weight)
      }
      
      readout_values <- c(readout_values, weighted_sum)
    } else {
      warning(sprintf("Readout %s not found in attractors", readout))
      readout_values <- c(readout_values, NA)
    }
  }
  
  names(readout_values) <- readouts
  return(readout_values)
}

#' Format data as proper MIDAS file with validation
#' @param midas_data Raw MIDAS data frame
#' @param stimuli Vector of stimulus node names
#' @param inhibitors Vector of inhibitor node names  
#' @param readouts Vector of readout node names
#' @return Properly formatted MIDAS data frame
format_midas <- function(midas_data, stimuli, inhibitors, readouts) {
  
  # Validate inputs
  if (nrow(midas_data) == 0) {
    stop("No MIDAS data to format")
  }
  
  # Create proper column names
  id_col <- "ID:treatment"
  stim_cols <- paste0("TR:", stimuli)
  inhib_cols <- paste0("TR:", inhibitors)
  readout_cols <- paste0("DA:", readouts)
  
  # Ensure all required columns exist
  all_cols <- c("ID", stimuli, inhibitors, readouts)
  missing_cols <- setdiff(all_cols, colnames(midas_data))
  
  if (length(missing_cols) > 0) {
    warning(sprintf("Missing columns in MIDAS data: %s", paste(missing_cols, collapse = ", ")))
    # Add missing columns with NA values
    for (col in missing_cols) {
      midas_data[[col]] <- NA
    }
  }
  
  # Reorder and rename columns
  final_midas <- midas_data[, c("ID", stimuli, inhibitors, readouts), drop = FALSE]
  colnames(final_midas) <- c(id_col, stim_cols, inhib_cols, readout_cols)
  
  # Convert to numeric where appropriate (except ID column)
  for (i in 2:ncol(final_midas)) {
    final_midas[, i] <- as.numeric(final_midas[, i])
  }
  
  return(final_midas)
}

#' Save MIDAS file with metadata
#' @param midas_data Formatted MIDAS data frame
#' @param filename Output filename
#' @param add_metadata Logical, whether to add metadata as comments
save_midas <- function(midas_data, filename, add_metadata = TRUE) {
  
  if (add_metadata) {
    # Create metadata header
    metadata <- c(
      "# MIDAS format data generated by BoolNet simulation",
      paste("# Generated on:", Sys.time()),
      paste("# Number of experiments:", nrow(midas_data)),
      paste("# Number of perturbations:", sum(grepl("^TR:", colnames(midas_data)))),
      paste("# Number of readouts:", sum(grepl("^DA:", colnames(midas_data)))),
      "#"
    )
    
    # Write metadata and data
    writeLines(metadata, filename)
    write.table(midas_data, filename, append = TRUE, sep = ",", 
                row.names = FALSE, quote = FALSE)
  } else {
    write.csv(midas_data, filename, row.names = FALSE, quote = FALSE)
  }
  
  cat(sprintf("MIDAS data saved to: %s\n", filename))
}

# Command-line interface
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-d", "--dataset"), type="character", default="tcell",
              help="Dataset name or path to model file", metavar="STRING"),
  make_option(c("-m", "--method"), type="character", default="random",
              help="Simulation method: random or exhaustive [default %default]", metavar="STRING"),
  make_option(c("-n", "--nsim"), type="integer", default=1000,
              help="Number of simulation iterations [default %default]", metavar="INT"),
  make_option(c("-c", "--combinations"), action="store_true", default=FALSE,
              help="Include combination experiments"),
  make_option(c("-s", "--seed"), type="integer", default=44,
              help="Random seed [default %default]", metavar="INT")
)

parser <- OptionParser(option_list=option_list,
                       description = "Generate MIDAS data from Boolean network models")

# Only parse arguments if script is run from command line
if (sys.nframe() == 0) {
  opt <- parse_args(parser)
  
  message("=== MIDAS Data Generation Pipeline ===")
  message(sprintf("Dataset: %s", opt$dataset))
  message(sprintf("Method: %s", opt$method))
  message(sprintf("Simulations: %d", opt$nsim))
  message(sprintf("Include combinations: %s", opt$combinations))
  message(sprintf("Random seed: %d", opt$seed))
  
  set.seed(opt$seed)
  
  # Load network based on dataset
  if (opt$dataset == "tcell") {
    sbml_file <- "data/T-Cell/T-Cell.sbml"
    bnet_file <- sub(".sbml$", ".bnet", sbml_file)
    midas_file <- "data/T-Cell/T-Cell_simulated.csv"
    
    if (!file.exists(sbml_file)) {
      stop(sprintf("SBML file not found: %s", sbml_file))
    }
    
    # Convert and load model
    convertSBML(sbml_file)
    model <- loadSBML(sbml_file)
    
  } else if (file.exists(opt$dataset)) {
    # Load from file
    if (grepl("\\.sbml$", opt$dataset)) {
      convertSBML(opt$dataset)
      model <- loadSBML(opt$dataset)
    } else if (grepl("\\.(rds|RDS)$", opt$dataset)) {
      model <- readRDS(opt$dataset)
    } else if (grepl("\\.(rdata|RData)$", opt$dataset)) {
      load(opt$dataset)  # Assumes model object is named 'model'
    } else {
      stop("Unsupported file format. Use .sbml, .rds, or .RData")
    }
    
    midas_file <- paste0(tools::file_path_sans_ext(opt$dataset), "_simulated.csv")
  } else {
    stop(sprintf("Dataset not recognized: %s", opt$dataset))
  }
  
  # Identify network components
  message("\n=== Analyzing Network Structure ===")
  node_types <- identify_node_types(model)
  
  message(sprintf("External nodes: %s", paste(node_types$external, collapse = ", ")))
  message(sprintf("Stimuli (%d): %s", length(node_types$stimuli), paste(node_types$stimuli, collapse = ", ")))
  message(sprintf("Inhibitors (%d): %s", length(node_types$inhibitors), paste(node_types$inhibitors, collapse = ", ")))
  message(sprintf("Readouts (%d): %s", length(node_types$readouts), paste(node_types$readouts, collapse = ", ")))
  
  # Design experiments
  message("\n=== Designing Experiments ===")
  experiments <- design_experiments(node_types$stimuli, node_types$inhibitors, 
                                   node_types$readouts, opt$combinations)
  message(sprintf("Total experiments: %d", length(experiments)))
  
  # Generate MIDAS data
  message("\n=== Running Simulations ===")
  midas_raw <- generate_midas_data(model, experiments, opt$nsim, opt$method)
  
  # Format and save
  message("\n=== Formatting and Saving Results ===")
  midas_formatted <- format_midas(midas_raw, node_types$stimuli, 
                                 node_types$inhibitors, node_types$readouts)
  save_midas(midas_formatted, midas_file)
  
  message(sprintf("\n=== Pipeline Complete ==="))
  message(sprintf("Results saved to: %s", midas_file))
  message(sprintf("Generated %d experiments with %d perturbations and %d readouts", 
                  nrow(midas_formatted),
                  length(node_types$stimuli) + length(node_types$inhibitors),
                  length(node_types$readouts)))
}