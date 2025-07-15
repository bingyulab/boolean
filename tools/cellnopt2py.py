#!/usr/bin/env python3
"""
CellNOpt Example using RPy2
This script demonstrates how to use CellNOpt with RPy2 for Boolean network optimization.
"""

import os
import sys
from rpy2.robjects.packages import importr
from rpy2.robjects import r


dataset_map = {
    "toy": ("ToyModel", "ToyModel.sif", "ToyModel.RData", "ToyModel.csv", "ToyModel.bnet"),
    "apoptosis": ("Apoptosis", "Apoptosis.sif", "Apoptosis.RData", "Apoptosis.csv", "Apoptosis.bnet"),
    "dream": ("DREAMmodel", "DreamModel.sif", "DreamModel.RData", "DreamModel.csv", "DreamModel.bnet"),
    "T-Cell": ("T-Cell", "TCell.sif", "TCell.RData", "TCell.csv", "TCell.bnet"),
}
        
class CellNOptAnalyzer:
    """
    A Python class to interface with CellNOpt R package using RPy2
    """

    def __init__(self, dataset="toy", ChangePct=0.1):
        """Initialize the CellNOpt analyzer"""
        self.cellnoptr = importr('CellNOptR')
        self.PERTUB = True if ChangePct > 0 else False
        self.ChangePct = ChangePct
        self.parse(dataset)
        print("CellNOptR loaded successfully")
    
    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name = dataset_map[dataset]
        self.dataset = base_name
        self.filename = "0_Modified" if not self.PERTUB else f"{self.ChangePct * 100:.0f}_Modified"
        self.input_path = os.path.join("data", base_name, self.filename)
        self.sif_file = os.path.join(self.input_path, sif_name)
        self.data_file = os.path.join(self.input_path, rdata_name)
        self.midas_file = os.path.join(self.input_path, midas_name)
        self.output_file = os.path.join("output/cellnopt", self.dataset, self.filename)       
        
        self.GD_MODEL = os.path.join("data", base_name, bnet_name)
        
        if not os.path.exists(self.output_file):
            os.makedirs(self.output_file, exist_ok=True)
    
    def load_network(self):
        """
        Load a network from file
        Args:
            file (str): Path to file
        Returns:
            R object containing the network
        """
        print(f"Loading network from {self.sif_file}")

        r(f'pknmodel <- readSIF("{self.sif_file}")')
        return r('pknmodel')
    
    def load_data(self):
        """
        Load experimental data from MIDAS file
        Args:
            midas_file (str): Path to MIDAS file
        Returns:
            R object containing the experimental data
        """        
        output_file = os.path.join(self.output_file, f"{self.dataset}_CNOlist.pdf")
        print(f"Loading data from {self.midas_file}, and saving to {output_file}")
        r(f'cnolist <- makeCNOlist(readMIDAS("{self.midas_file}", verbose=TRUE), subfield=FALSE)')
        r(f'plotCNOlistPDF(CNOlist=CNOlist(cnolist), filename="{output_file}")') 
        return r('cnolist')
    
    def plot_network(self, cnolistPB=None, output_file="output/cellnopt/network_plot"):
        """
        Plot the network and save to PDF
        Args:
            output_file (str): Path to output PDF file
        """
        from rpy2.robjects import r
        from rpy2.robjects.packages import importr
        output_file = os.path.join(self.output_file, f"{self.dataset}_network_plot")
        print(f"Plotting network and saving to {output_file}")
        r(f'library(CellNOptR)')
        r(f'library(CNORdt)')
        r(f'library(Rgraphviz)')
        r(f'library(RBGL)')
        r(f'model <- readSBMLQual(file)')
        if cnolistPB is not None:
            r(f'plotModel(model, filename="{output_file}", output="SVG")')
        else:
            r(f'plotModel(model, cnolist, filename="{output_file}", output="SVG")')

    def preprocess_network(self):
        """
        Preprocess the network and data
        Finding and cutting the non observable and non controllable species
        Compressing the model
        Expanding the gates
        Returns:
            R object containing the preprocessed model
        """
        print("Preprocessing network...")
        r('model <- preprocessing(data = cnolist, model = pknmodel)')
        return r('model')

    def optimize_network(self, output_file, PLOT=True, method='ga', numSol=3, relGap=0.05):
        """
        Run Boolean network optimization
        Args:
            output_file (str): Directory to save output files within method directory, different for self.output_file
            PLOT (bool): Whether to plot the results
            method (str): Optimization method to use ('ga' or 'ilp')
        Returns:
            R object containing optimization results
        """
        print(f"Running optimization (method = {method}, number of solution ={numSol}, relGap={relGap})...")
        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')   
        if method == 'ga':  
            r(f'''
                resEcnolist <- residualError(cnolist);
                initBstring <- rep(1,length(model$reacID));
                maxGens=1000;
                stallGenMax=Inf;
                popSize=200;
                elitism=popSize/10;
                sizeFac=1e-04; 
                NAFac=1; 
                maxTime=Inf; 
                numStarts = 5;
                t <- system.time(opt_results<-gaBinaryT1(
                    CNOlist=cnolist,
                    model=model,
                    sizeFac=sizeFac,
                    initBstring=initBstring,
                    maxGens=maxGens, popSize=popSize,  elitism=elitism, 
                    stallGenMax=stallGenMax,
                    maxTime=maxTime, verbose=TRUE)
                )
                cat("Time taken for optimization:", t[3], "seconds\\n");
                optModel <- cutModel(model, opt_results$bString);
                simResults <- simulate_CNO(model=model,#_orig,
                            CNOlist=cnolist,
                            bString=opt_results$bString)                
            ''')
            if PLOT:
                r(f'''
                    cutAndPlot(model=model, bStrings=list(opt_results$bString),
                        CNOlist=cnolist,plotPDF=TRUE);
                    pdf("{output_file}/{self.dataset}_evolFitT1.pdf");
                    plotFit(optRes=opt_results);
                    dev.off();
                    plotModel(model, cnolist, bString=opt_results$bString, output="SVG", filename="{output_file}/{self.dataset}_mapback_evolFitT1_1.svg");
                    save(simResults,file=paste("{output_file}/{self.dataset}_evolSimRes.RData",sep=""))                    
                ''')
        elif method == 'ilp':   
            cplexPath = "/home/users/bjiang/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex" 
            r(f'''
                print("Running ILP optimization...");
                cnolist <- CNOlist(cnolist);
                resILPAll <- list();
                exclusionList <- NULL;
                cnolistReal <- cnolist;       
                writeMIDAS(CNOlist = cnolist, filename = "tempMD.csv", 
                        timeIndices = c(1, 2), overwrite = TRUE);
                md <- readMIDAS(MIDASfile = "tempMD.csv");
                file.remove("tempMD.csv");
                cnolist <- compatCNOlist(object = cnolist); 
                options(scipen = 10); 
                cnolist <- makeCNOlist(md,FALSE);                           
            ''')
            r(f'''  
                startILP <- Sys.time();
                print("Creating LP file and running ILP...");   
                resILP <- CellNOptR:::createAndRunILP(
                    model, md, cnolistReal, accountForModelSize = TRUE, 
                    sizeFac = 0.0001, mipGap=0, relGap={relGap}, 
                    timelimit=3600, cplexPath = "{cplexPath}", 
                    method = "quadratic", numSolutions = {numSol}, 
                    limitPop = 500, poolIntensity = 0, poolReplace = 2);
                endILP <- Sys.time();                
                CellNOptR:::cleanupILP();
                opt_results <- resILP;                
                optModel <- cutModel(model, opt_results$bitstringILP[[1]]);
                simResults <- simulate_CNO(model=model,#_orig,
                            CNOlist=cnolist,
                            bString=opt_results$bitstringILP[[1]])
            ''')
            if PLOT:
                r(f'''
                    cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(opt_results$bitstringILP[[1]]), plotPDF=TRUE)
                    plotModel(model, cnolist, bString=opt_results$bitstringILP[[1]], output="SVG", filename="{output_file}/{self.dataset}_ilpFitT1.svg");
                    save(simResults,file=paste("{output_file}/{self.dataset}_ilpSimRes.RData",sep=""))  
                ''')
        return r('optModel')          

    def save_results(self, output_dir):
        """
        Save optimization results
        Args:
            output_dir (str): Output directory
        """
        print(f"Saving results to {output_dir}/")        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')
        r(f'''
            if (ncol(optModel$notMat) == 0) {{
                n <- nrow(optModel$interMat)
                m <- ncol(optModel$interMat)
                optModel$notMat <- matrix(0, nrow=n, ncol=m)
                rownames(optModel$notMat) <- rownames(optModel$interMat)
                colnames(optModel$notMat) <- colnames(optModel$interMat)
            }}
        ''')
        sif_fname = os.path.join(output_dir, f"OPT_{self.dataset}.sif")
        rdata_fname = os.path.join(output_dir, f"OPT_{self.dataset}.RData")
        boolnet_fname = os.path.join(output_dir, f"OPT_{self.dataset}.bnet")
        
        print(f"Saving SIF to {sif_fname}...")
        
        r(f'''
            writeSIF(optModel, file="{sif_fname}", overwrite=TRUE)
        ''')
        print(f"Saving RData to {rdata_fname}...")
        r(f'''
            save(optModel, file="{rdata_fname}")
        ''')
        print(f"Saving BoolNet to {boolnet_fname}...")
        r(f'''            
            SIFToBoolNet(sifFile     = "{sif_fname}",
                        boolnetFile = "{boolnet_fname}",
                        CNOlist     = cnolist,
                        model       = optModel,
                        fixInputs   = FALSE,
                        preprocess  = TRUE,
                        ignoreAnds  = TRUE)
        ''')
        r(f'''
            namesFiles<-list(
                dataPlot="{output_dir}/ModelGraph.pdf",
                evolFitT1="{output_dir}/evolFitT1.pdf",
                evolFitT2=NA,
                simResultsT1="{output_dir}/SimResultsT1_1.pdf",
                simResultsT2=NA,
                scaffold="{output_dir}/Scaffold.sif",
                scaffoldDot="{output_dir}/Scaffold.dot",
                tscaffold="{output_dir}/TimesScaffold.EA",
                wscaffold="{output_dir}/weightsScaffold.EA",
                PKN="{output_dir}/PKN.sif",
                PKNdot="{output_dir}/PKN.dot",
                wPKN="{output_dir}/TimesPKN.EA",
                nPKN="{output_dir}/nodesPKN.NA")
        ''')
        r(f'''
            writeReport(
                modelOriginal=pknmodel,
                modelOpt=model,
                optimResT1=opt_results,
                optimResT2=NA,
                CNOlist=cnolist,
                directory="test",
                namesFiles=namesFiles,
                namesData=list(CNOlist="cnolist",model="Model"))
        ''')
        
    def evaluate_model(self, output_dir):
        """
        Evaluate the model using the CellNOptR package
        """
        print("Evaluating model...")
        r('library(here)')
        r('source(here::here("tools", "comparison.R"))')
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(output_dir, f"OPT_{self.dataset}.bnet")
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")
        r(f'''
          res <- compareNetwork(origFile = "{self.GD_MODEL}", modifiedFiles = list("{opt_fname}"))
          ''')
        res = r['res']
        from functions import parse_rpy2_results, analyze_attractor_performance
        results = parse_rpy2_results(res)

        analyze_attractor_performance(results[f"OPT_{self.dataset}.bnet"])
        # Now result_dict is a Python dictionary
        
        if 'ilp' in output_dir:
            results[fname]['bitstring'] = [r('opt_results$bitstringILP[[1]]')] * len(results[fname])
            results[fname]['total_time'] = [r('opt_results$cpu_time')] * len(results[fname])
        else:
            results[fname]['bitstring'] = [r('opt_results$bString')] * len(results[fname])
            results[fname]['total_time'] = [r('opt_results$total_time')] * len(results[fname])

        pd.DataFrame(results).to_csv(os.path.join(output_dir, "results.csv"), index=False)

        return results  
    
    def run_full_analysis(self, method="ga", numSol=3, relGap=0.05,):
        """
        Run complete CellNOpt analysis pipeline
        Args:
            file (str): Path to network file
            midas_file (str): Path to MIDAS data file
            output_dir (str): Output directory for results
            numSol (int): Number of solutions for ILP optimization
            relGap (float): Relative gap for ILP optimization
        """
        # Load network and data
        model = self.load_network()
        cnolist = self.load_data()
        
        # Preprocess
        self.preprocess_network()

        # Optimize
        output_file = os.path.join(self.output_file, method)                
        os.makedirs(output_file, exist_ok=True)
        print(f"Optimizing with method: {method}")
        opt_results = self.optimize_network(PLOT=True, method=method,
                                numSol=3, relGap=0.05, output_file=output_file)

        # Save results
        self.save_results(output_file)
        print(f"Results saved to {output_file}/")
        
        # Evaluate and cross-validate
        print("Evaluating and cross-validating the model...")
        results = self.evaluate_model(output_file)        
        print("Analysis completed successfully!")
        return model, cnolist, results

def main(dataset="toy", PERTUB=False, ChangePct=0.1, method="ga"):
    """Example usage of CellNOptAnalyzer"""
    
    # file = "data/apoptosis.xml"
    # For MIDAS file, you'll need to create or provide one
    # This is just an example - you'll need actual experimental data
    # midas_file = "data/experimental_data.csv"  # You need to create this

    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer(
        dataset=dataset,
        ChangePct=ChangePct
    )

    model, cnolist, results = analyzer.run_full_analysis(
        method=method,
        numSol=3,
        relGap=0.05
    )
    return model, cnolist, results

if __name__ == "__main__":
    # model, cnolist, results = main(Test=True, method="ilp")
    model, cnolist, results = main(dataset="toy", ChangePct=0.1, method="ga")