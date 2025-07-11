#!/usr/bin/env python3
"""
CellNOpt Example using RPy2
This script demonstrates how to use CellNOpt with RPy2 for Boolean network optimization.
"""

import os
import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


class CellNOptAnalyzer:
    """
    A Python class to interface with CellNOpt R package using RPy2
    """
    
    def __init__(self, file, midas_file=None):
        """Initialize the CellNOpt analyzer"""
        self.cellnoptr = importr('CellNOptR')
        self.file = file
        self.midas_file = midas_file
        if file is None:
            self.filename = "OriginalModel"
        else:
            self.filename = os.path.splitext(os.path.basename(file))[0]
        self.output_file = "output/cellnopt"
        if not os.path.exists(self.output_file):
            os.makedirs(self.output_file, exist_ok=True)
        print("CellNOptR loaded successfully")
    
    def load_network(self):
        """
        Load a network from file
        Args:
            file (str): Path to file
        Returns:
            R object containing the network
        """
        print(f"Loading network from {self.file}")
        if self.file is None:
            robjects.r(f'data(ToyModel, package="CellNOptR")')
            robjects.r('pknmodel <- ToyModel')
        elif self.file.endswith('.sif'):
            robjects.r(f'pknmodel <- readSIF("{self.file}")')
        elif self.file.endswith('.xml'):
            robjects.r(f'pknmodel <- readSBMLQual("{self.file}")')
        elif self.file.endswith('.RData'):
            robjects.r(f'load("{self.file}")')
            robjects.r('pknmodel <- model')

        return robjects.r('pknmodel')
    
    def load_data(self):
        """
        Load experimental data from MIDAS file
        Args:
            midas_file (str): Path to MIDAS file
        Returns:
            R object containing the experimental data
        """        
        output_file = os.path.join(self.output_file, f"{self.filename}_CNOlist.pdf")
        print(f"Loading data from {self.midas_file}, and saving to {output_file}")
        if self.midas_file is None:
            robjects.r(f'data(CNOlistToy, package="CellNOptR")')
            robjects.r('cnolist <- CNOlistToy')  
            robjects.r(f'plotCNOlistPDF(CNOlist=cnolist, filename="{output_file}")')  # Example timepoints
        else:  
            robjects.r(f'data <- readMIDAS("{midas_file}")')
            robjects.r(f'cnolist <- makecnolist(data, subfield=FALSE, verbose=FALSE)')  # Remove NA timepoints
            robjects.r('cnolist <- CNOlist(cnolist)')  
            robjects.r(f'plotCNOlistPDF(CNOlist=cnolist, filename="{output_file}")')  # Example timepoints
        return robjects.r('cnolist')
    
    def plot_network(self, cnolistPB=None, output_file="output/cellnopt/network_plot"):
        """
        Plot the network and save to PDF
        Args:
            output_file (str): Path to output PDF file
        """
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        output_file = os.path.join(self.output_file, f"{self.filename}_network_plot")
        print(f"Plotting network and saving to {output_file}")
        robjects.r(f'library(CellNOptR)')
        robjects.r(f'library(CNORdt)')
        robjects.r(f'library(Rgraphviz)')
        robjects.r(f'library(RBGL)')
        robjects.r(f'model <- readSBMLQual(file)')
        if cnolistPB is not None:
            robjects.r(f'plotModel(model, filename="{output_file}", output="SVG")')
        else:
            robjects.r(f'plotModel(model, cnolist, filename="{output_file}", output="SVG")')

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
        robjects.r('model <- preprocessing(data = cnolist, model = pknmodel)')
        return robjects.r('model')

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
            robjects.r(f'''
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
                save(res,file=paste("{output_file}/{self.filename}_evolSimRes.RData",sep=""))                    

            ''')
            if PLOT:
                robjects.r(f'''
                    cutAndPlot(model=model, bStrings=list(opt_results$bString),
                        CNOlist=cnolist,plotPDF=TRUE);
                    pdf("{output_file}/{self.filename}_evolFitT1.pdf");
                    plotFit(optRes=opt_results);
                    dev.off();
                    plotModel(model, cnolist, bString=opt_results$bString);
                    bs <- mapBack(model, pknmodel, opt_results$bString);
                    plotModel(pknmodel, cnolist, bs, compressed=model$speciesCompressed, output="SVG", filename="{output_file}/{self.filename}_evolFitT1_{i}.svg") ;
                ''')
        elif method == 'ilp':   
            cplexPath = "/home/users/bjiang/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex" 
            robjects.r(f'''
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
            robjects.r(f'''  
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
                robjects.r(f'''
                    cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(opt_results$bitstringILP[[1]]), plotPDF=TRUE)
                    plotModel(model, cnolist, bString=opt_results$bitstringILP[[1]], output="SVG", filename="{output_file}/{self.filename}_ilpFitT1_1.svg");
                    save(res,file=paste("{output_file}/{self.filename}_ilpSimRes.RData",sep=""))  
                ''')
        return robjects.r('optModel')      
    
    def evaluate_model(self, method='ga'):
        """
        Evaluate the model using the CellNOptR package
        """
        print("Evaluating model...")
        r('library(here)')
        r('source(here::here("tools", "comparison.R"))')      
        
        ori_fname = os.path.join(output_dir, f"{self.filename}.txt")
        opt_fname = os.path.join(output_dir, f"OPT_{self.filename}.txt")
        r(f'''
          res <- compareNetworks(origFile = {ori_fname}, modifiedFiles = list({opt_fname}))
          ''')
        res = r['res']

        from tools.functions import rlist_to_pydict
        result_dict = rlist_to_pydict(res)

        # Now result_dict is a Python dictionary
        return result_dict        
    
    def prepare_cnodata(self):
        # Select timepoints and cut cnolist
        r('selectedTime <- c(0,10)')
        r('cnodata_prep <- cutcnolist(cnodata, model = model, cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))')
        r('plot(cnodata_prep)')
        
    def cross_validate(self):
        # Register parallel backend
        r('doParallel::registerDoParallel(cores=3)')
        # 10-fold cross-validation with different types and parallelization
        print("Cross-validation (datapoint, parallel=TRUE):")
        r('system.time({R1 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="datapoint", nfolds=10, parallel = TRUE)})')
        print("Cross-validation (experiment, parallel=TRUE):")
        r('system.time({R2 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="experiment", nfolds=10, parallel = TRUE)})')
        print("Cross-validation (observable, parallel=TRUE):")
        r('system.time({R3 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="observable", nfolds=10, parallel = TRUE)})')
        print("Cross-validation (datapoint, parallel=FALSE):")
        r('system.time({R4 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="datapoint", nfolds=10, parallel = FALSE)})')
        print("Cross-validation (experiment, parallel=FALSE):")
        r('system.time({R5 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="experiment", nfolds=10, parallel = FALSE)})')
        print("Cross-validation (observable, parallel=FALSE):")
        r('system.time({R6 <- crossvalidateBoolean(CNOlist = cnodata_prep, model = model, type="observable", nfolds=10, parallel = FALSE)})')
        print("Done.")

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
        sif_fname = os.path.join(output_dir, f"OPT_{self.filename}.sif")
        rdata_fname = os.path.join(output_dir, f"OPT_{self.filename}.RData")
        boolnet_fname = os.path.join(output_dir, f"OPT_{self.filename}.txt")
        
        # filename <- tools::file_path_sans_ext(basename(file))
        sif_fname    <- file.path(output_file,     sprintf("OPT_%s.sif",  self.filename))
        rdata_fname  <- file.path(output_file,     sprintf("OPT_%s.RData",self.filename))
        boolnet_fname<- file.path(out_dir_boolnet, sprintf("OPT_%s.txt",  self.filename))
  
        robjects.r(f'''
            save(optModel, file={rdata_fname})
            writeSIF(optModel, file={sif_fname}, overwrite=TRUE)
            message("Wrote:\n - RData → ", {rdata_fname}, "\n - SIF   → ", {sif_fname}, "\n")
        ''')
        robjects.r(f'''            
            SIFToBoolNet(sifFile     = {sif_fname},
                        boolnetFile = {boolnet_fname},
                        CNOlist     = CNOlist,
                        model       = optModel,
                        fixInputs   = TRUE,
                        preprocess  = TRUE,
                        ignoreAnds  = TRUE)
            message("BoolNet file written to: ", {boolnet_fname})
        ''')
        robjects.r(f'''
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
        robjects.r(f'''
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

        # Evaluate and cross-validate
        print("Evaluating and cross-validating the model...")
        self.evaluate_model()
        
        # Save results
        self.save_results(output_file)
        
        print("Analysis completed successfully!")
        return model, cnolist, opt_results

def main(file, midas_file=None, method="ga"):
    """Example usage of CellNOptAnalyzer"""
    
    # file = "data/apoptosis.xml"
    # For MIDAS file, you'll need to create or provide one
    # This is just an example - you'll need actual experimental data
    # midas_file = "data/experimental_data.csv"  # You need to create this

    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer(
        file=file,
        midas_file=midas_file)

    model, cnolist, opt_results = analyzer.run_full_analysis(
        method=method,
        numSol=3,
        relGap=0.05
    )
    return model, cnolist, opt_results

if __name__ == "__main__":
    # model, cnolist, opt_results = main(Test=True, method="ilp")
    model, cnolist, opt_results = main(
        "./output/ModifiedToyModel_10.RData", 
        midas_file=None, method="ilp")