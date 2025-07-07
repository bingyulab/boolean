#!/usr/bin/env python3
"""
CellNOpt Example using RPy2
This script demonstrates how to use CellNOpt with RPy2 for Boolean network optimization.
"""

import os
import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

# Activate automatic numpy to R conversion
rpy2.robjects.numpy2ri.activate()

class CellNOptAnalyzer:
    """
    A Python class to interface with CellNOpt R package using RPy2
    """
    
    def __init__(self):
        """Initialize the CellNOpt analyzer"""
        try:
            # Import required R packages
            self.base = importr('base')
            self.cellnoptr = importr('CellNOptR')
            print("CellNOptR loaded successfully")
        except Exception as e:
            print(f"Error loading R packages: {e}")
            print("Make sure you have installed CellNOptR in R:")
            print("BiocManager::install('CellNOptR')")
            sys.exit(1)
    
    def load_network(self, sif_file):
        """
        Load a network from SIF file
        Args:
            sif_file (str): Path to SIF file
        Returns:
            R object containing the network
        """
        if not os.path.exists(sif_file):
            raise FileNotFoundError(f"SIF file not found: {sif_file}")
        
        print(f"Loading network from {sif_file}")
        robjects.r(f'pknmodel <- readSIF("{sif_file}")')
        return robjects.r('pknmodel')
    
    def load_data(self, midas_file):
        """
        Load experimental data from MIDAS file
        Args:
            midas_file (str): Path to MIDAS file
        Returns:
            R object containing the experimental data
        """
        if not os.path.exists(midas_file):
            raise FileNotFoundError(f"MIDAS file not found: {midas_file}")
        
        print(f"Loading data from {midas_file}")
        robjects.r(f'CNOlist <- readMIDAS("{midas_file}")')
        return robjects.r('CNOlist')
    
    def plot_network(self, file, CNOlistPB=None, output_file="output/network_plot.pdf"):
        """
        Plot the network and save to PDF
        Args:
            output_file (str): Path to output PDF file
        """
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        robjects.r(f'library(CellNOptR)')
        robjects.r(f'library(CNORdt)')
        robjects.r(f'library(Rgraphviz)')
        robjects.r(f'library(RBGL)')
        robjects.r(f'model <- readSBMLQual("{file}")')
        robjects.r(f'data(CNOlistPB, package="CNORdt")')
        if CNOlistPB is not None:
            robjects.r(f'CNOlistPB <- readMIDAS("{CNOlistPB}")')
        else:
            robjects.r('CNOlistPB <- CNOlist')
            robjects.r(f'plotModel(model, CNOlistPB, filename="{output_file}", output="SVG")')
    
    def preprocess_network(self):
        """
        Preprocess the network and data
        Returns:
            R object containing the preprocessed model
        """
        print("Preprocessing network...")
        robjects.r('model <- preprocessing(CNOlist, pknmodel)')
        return robjects.r('model')
    
    def optimize_boolean_network(self, max_generations=500, population_size=50):
        """
        Run Boolean network optimization
        Args:
            max_generations (int): Maximum number of generations for GA
            population_size (int): Population size for GA
        Returns:
            R object containing optimization results
        """
        print(f"Running optimization (maxGens={max_generations}, popSize={population_size})...")
        robjects.r(f'''
        opt_results <- gaBinaryT1(
            CNOlist = CNOlist, 
            model = model, 
            maxGens = {max_generations}, 
            popSize = {population_size},
            verbose = TRUE
        )
        ''')
        return robjects.r('opt_results')
    
    def save_results(self, output_dir="output"):
        """
        Save optimization results
        Args:
            output_dir (str): Output directory
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print(f"Saving results to {output_dir}/")
        robjects.r(f'''
        # Save optimized network
        writeNetwork(model=model, optimResT1=opt_results, 
                    filename="{output_dir}/optimized_network")
        
        # Save results summary
        capture.output(print(opt_results), file="{output_dir}/optimization_summary.txt")
        ''')
    
    def run_full_analysis(self, sif_file, midas_file, output_dir="output", 
                         max_generations=500, population_size=50):
        """
        Run complete CellNOpt analysis pipeline
        Args:
            sif_file (str): Path to SIF network file
            midas_file (str): Path to MIDAS data file
            output_dir (str): Output directory for results
            max_generations (int): Maximum generations for optimization
            population_size (int): Population size for optimization
        """
        try:
            # Load network and data
            self.load_network(sif_file)
            self.load_data(midas_file)
            
            # Preprocess
            self.preprocess_network()
            
            # Optimize
            self.optimize_boolean_network(max_generations, population_size)
            
            # Save results
            self.save_results(output_dir)
            
            print("Analysis completed successfully!")
            
        except Exception as e:
            print(f"Error during analysis: {e}")
            raise

def main():
    """Example usage of CellNOptAnalyzer"""
    
    # Check if we have the required files
    sif_file = "data/apoptosis.sif"
    
    # For MIDAS file, you'll need to create or provide one
    # This is just an example - you'll need actual experimental data
    midas_file = "data/experimental_data.csv"  # You need to create this
    
    if not os.path.exists(sif_file):
        print(f"Example SIF file not found: {sif_file}")
        print("Please make sure you have a network file in SIF format")
        return
    
    if not os.path.exists(midas_file):
        print(f"MIDAS data file not found: {midas_file}")
        print("You need to create experimental data in MIDAS format")
        print("See CellNOpt documentation for format details")
        return
    
    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer()
    
    try:
        analyzer.run_full_analysis(
            sif_file=sif_file,
            midas_file=midas_file,
            output_dir="output/cellnopt_results",
            max_generations=100,  # Reduced for testing
            population_size=20    # Reduced for testing
        )
    except Exception as e:
        print(f"Analysis failed: {e}")

if __name__ == "__main__":
    main()

