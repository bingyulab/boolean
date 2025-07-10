#!/usr/bin/env Rscript
# Install CellNOpt and related packages

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

cat("Installing MEIGOR packages...\n") # browseVignettes("MEIGOR")
BiocManager::install("MEIGOR")

# Install Bioconductor packages
cat("Installing CellNOpt packages...\n")

packages_to_install <- c(
    "CellNOptR",      # Main CellNOpt package
    "CNORdt",         # Discrete time modeling
    "CNORfeeder",     # Network expansion
    "CNORfuzzy",      # Fuzzy logic modeling
    "graph",          # Graph utilities
    "Rgraphviz",       # Graph visualization (if needed)
    "Cairo",
    "BoolNet",
    "here"
)

for (pkg in packages_to_install) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        cat(sprintf("✓ %s installed successfully\n", pkg))
    }, error = function(e) {
        cat(sprintf("✗ Failed to install %s: %s\n", pkg, e$message))
    })
}

# Test installation
cat("Testing CellNOptR installation...\n")
tryCatch({
    library(CellNOptR)
    cat("✓ CellNOptR loaded successfully\n")
    
    # Print version info
    cat(sprintf("CellNOptR version: %s\n", packageVersion("CellNOptR")))
    
}, error = function(e) {
    cat(sprintf("✗ Failed to load CellNOptR: %s\n", e$message))
})

cat("Installation complete!\n")
