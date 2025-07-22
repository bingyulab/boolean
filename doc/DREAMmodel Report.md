```r
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
summarizeDreamData <- function(model, CNOl) {
  # 1) Network (BoolNet) summary
  n_nodes       <- length(model$namesSpecies)
  n_reactions   <- length(model$reacID)
  n_activations <- sum(model$interMat == 1)
  n_inhibitions <- sum(model$notMat   == 1)
  
  cat("=== Network summary ===\n")
  cat("• Nodes (species)         :", n_nodes, "\n")
  cat("• Reactions (edges)       :", n_reactions, "\n")
  cat("  – Activating edges      :", n_activations, "\n")
  cat("  – Inhibitory edges      :", n_inhibitions, "\n\n")
  
  # 2) Experimental data (CNOlist) summary
  n_cues       <- length(CNOl$namesCues)
  n_stimuli    <- length(CNOl$namesStimuli)
  n_inhibitors <- length(CNOl$namesInhibitors)
  n_signals    <- length(CNOl$namesSignals)
  n_times      <- length(CNOl$timeSignals)
  n_samples    <- nrow(CNOl$valueCues)
  
  cat("=== Experimental data summary ===\n")
  cat("• Cues (inputs measured)    :", n_cues, "\n")
  cat("• Stimuli applied           :", n_stimuli, "\n")
  cat("• Inhibitors applied        :", n_inhibitors, "\n")
  cat("• Signals read out          :", n_signals, "\n")
  cat("• Time‐points per experiment:", n_times, "\n")
  cat("• Total experiments (rows)  :", n_samples, "\n")
}
```
=== Network summary ===
• Nodes (species)         : 40 
• Reactions (edges)       : 58 
  – Activating edges      : 58 
  – Inhibitory edges      : 2 

=== Experimental data summary ===
• Cues (inputs measured)    : 8 
• Stimuli applied           : 4 
• Inhibitors applied        : 4 
• Signals read out          : 7 
• Time‐points per experiment: 2 
• Total experiments (rows)  : 25