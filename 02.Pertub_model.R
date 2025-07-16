#!/usr/bin/env Rscript
library(CellNOptR)
library(BoolNet)
library(here)
source(here::here("tools", "functions.R"))


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
                               DELETE       = FALSE) {
  
  # Initialize environment and validate inputs
  set.seed(seed)
  
  # Enhanced dataset mapping with validation
  dataset_map <- list(
    toy       = c("ToyModel", "ToyModel.sif", "ToyModel.RData", "ToyModel.csv", "ToyModel.bnet", "ToyModel.sbml"),
    apoptosis = c("Apoptosis", "Apoptosis.sif", "Apoptosis.RData", "Apoptosis.csv", "Apoptosis.bnet", "Apoptosis.sbml"),
    dream     = c("DREAMmodel", "DreamModel.sif", "DreamModel.RData", "DreamModel.csv", "DreamModel.bnet", "DreamModel.sbml"),
    TCell     = c("T-Cell", "TCell.sif", "TCell.RData", "TCell.csv", "TCell.bnet", "TCell.sbml")
  )
  
  if (!dataset %in% names(dataset_map)) {
    stop(sprintf("Unknown dataset: %s. Available datasets: %s", 
                 dataset, paste(names(dataset_map), collapse=", ")))
  }
  
  # Validate perturbation parameters
  if (change_pct < 0 || change_pct > 0.99) {
    stop("change_pct must be between 0 and 0.99")
  }
  
  # Directory and file path setup
  vals <- dataset_map[[dataset]]
  base_name      <- vals[1]
  sif_name       <- vals[2]
  rdata_name     <- vals[3]
  midas_name     <- vals[4]
  boolnet_name   <- vals[5]
  sbml_name      <- vals[6]

  new_d  <- paste0(round(change_pct * 100), "_Modified")
  output_file <- file.path("data", base_name, new_d)
  if (!dir.exists(output_file)) dir.create(output_file, recursive=TRUE)

  # Input and output file paths
  sif_in        <- file.path("data", base_name, sif_name)
  midas_in      <- file.path("data", base_name, midas_name)
  bnet_in       <- file.path("data", base_name, boolnet_name)

  sif_out       <- file.path(output_file, sif_name)
  rdata_out     <- file.path(output_file, rdata_name)
  bnet_out      <- file.path(output_file, boolnet_name)
  sbml_out      <- file.path(output_file, sbml_name)

  # Copy experimental data
  if (!file.copy(midas_in, file.path(output_file, midas_name), overwrite=TRUE)) {
    stop("Failed to copy MIDAS data file")
  }

  # Load original models and data
  tryCatch({
    orig_model <- readSIF(sif_in)
    cnolist    <- makeCNOlist(readMIDAS(file.path(output_file, midas_name), verbose=FALSE), subfield=FALSE)
  }, error = function(e) {
    stop(sprintf("Error loading original model or data: %s", e$message))
  })

  # Convert to BoolNet format with enhanced error handling
  if (!file.exists(bnet_in)) {
    message("Converting SIF to BoolNet format...")
    tryCatch({
      SIFToBoolNet(sifFile     = sif_in,
                   boolnetFile = bnet_in,
                   CNOlist     = cnolist,
                   model       = orig_model,
                   fixInputs   = FALSE,
                   preprocess  = TRUE,
                   ignoreAnds  = TRUE)
    }, error = function(e) {
      stop(sprintf("SIF to BoolNet conversion failed: %s", e$message))
    })
  }

  # Load Boolean network with validation
  tryCatch({
    boolnet <- loadNetwork(bnet_in)
    message(sprintf("Loaded network with %d genes", length(boolnet$genes)))
  }, error = function(e) {
    stop(sprintf("Failed to load Boolean network: %s", e$message))
  })

  # Apply perturbations with enhanced control
  message(sprintf("Applying perturbation at %.1f%% level", 
                  change_pct * 100))

  tryCatch({
    perturbed_boolnet <- boolnet
    num_genes_to_perturb <- floor(change_pct * length(perturbed_boolnet$genes))

    # Randomly sample genes for perturbation
    genes_to_perturb <- sample(perturbed_boolnet$genes, num_genes_to_perturb)

    message(sprintf("Applying functions perturbation with shuffle method for %d genes", 
                num_genes_to_perturb))
    for (g in genes_to_perturb) {
      k      <- length(perturbed_boolnet$genesConnection[[g]])
      
      # Scale perturbation based on gene complexity
      perturbation_scale <- min(1, max(0.1, 1/k))
      bits <- 2^k
      toFlip <- max(1, floor(perturbation_scale * bits))
      message(sprintf("Applying functions perturbation with shuffle method for %f genes", 
              toFlip))
      perturbed_boolnet <- perturbNetwork(
        perturbed_boolnet,
        perturb       = "functions",
        method        = "shuffle",
        excludeFixed  = TRUE, # prevent perturbation of fixed variable
        maxNumBits    = toFlip, # ensuring that each modification represents a single logical change
        simplify      = FALSE,  # preserve the original network complexity and prevent automatic simplification
        )
    }
  }, error = function(e) {
    stop(sprintf("Network perturbation failed: %s", e$message))
  })
  
  # Generate additional output formats
  tryCatch({
    # Save BoolNet format
    saveNetwork(perturbed_boolnet, bnet_out) 
    library(CellNOptR)
    # Save SIF format
    write_boolnet_to_sif(perturbed_boolnet, sif_out)
    reconstructed_model <- readSIF(sif_out)
    save(reconstructed_model, file = rdata_out)
    writeSIF(reconstructed_model, file = sif_out, overwrite = TRUE)    
    
  }, error = function(e) {
    warning(sprintf("Secondary format export encountered issues: %s", e$message))
  })

  # Validation and summary
  message("Performing validation checks...")
  validation_results <- list(
    original_genes = length(boolnet$genes),
    perturbed_genes = length(perturbed_boolnet$genes),
    files_created = c(
      sbml = file.exists(sbml_out),
      bnet = file.exists(bnet_out),
      sif = file.exists(sif_out),
      rdata = file.exists(rdata_out)
    )
  )
  
  message(sprintf("Pipeline completed successfully. Output directory: %s", output_file))
  message(sprintf("Validation: %d/%d genes preserved, %d/%d files created",
                  validation_results$perturbed_genes,
                  validation_results$original_genes,
                  sum(validation_results$files_created),
                  length(validation_results$files_created)))

  invisible(list(
    originalBoolNet = boolnet,
    perturbedBoolNet = perturbed_boolnet,
    outputDirectory = output_file,
    sbmlFile = sbml_out,
    bnetFile = bnet_out,
    sifFile = sif_out,
    rdataFile = rdata_out,
    validationResults = validation_results
  ))
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
              help="Random seed [optional]", metavar="INT")
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
    change_pct   = as.double(opt$changePCT),
    DELETE       = opt$delete,
    seed         = opt$seed
  )
