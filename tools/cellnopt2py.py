#!/usr/bin/env python3
"""
CellNOpt Example using RPy2
This script demonstrates how to use CellNOpt with RPy2 for Boolean network optimization.
"""

import os
import sys
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from tools.config import dataset_map
from tools.comparison import limit_float, AttractorAnalysis

        
class CellNOptAnalyzer:
    """
    A Python class to interface with CellNOpt R package using RPy2
    """

    def __init__(self, dataset="toy", manager=None):
        """Initialize the CellNOpt analyzer"""
        self.cellnoptr = importr('CellNOptR')
        self.PERTUB = True if manager.change_percent > 0 else False
        self.ChangePct = manager.change_percent
        self.ilp_config = manager.get_ilp_config()
        self.ga_config = manager.get_ga_config()
        self.parse(dataset)
        print("CellNOptR loaded successfully")
    
    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name, _ = dataset_map[dataset]
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

    def optimize_network(self, output_file, PLOT=True, method='ga'):
        """
        Run Boolean network optimization
        Args:
            output_file (str): Directory to save output files within method directory, different for self.output_file
            PLOT (bool): Whether to plot the results
            method (str): Optimization method to use ('ga' or 'ilp')
        Returns:
            R object containing optimization results
        """
        print(f"Running optimization (method = {method}...")
        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')   

        if method == 'ga':  
            r(f'''
                resEcnolist <- residualError(cnolist);
                initBstring <- rep(1,length(model$reacID));
                t <- system.time(opt_results<-gaBinaryT1(
                    CNOlist=cnolist, model=model, initBstring=initBstring,
                    maxGens={self.ga_config["maxGens"]}, sizeFac={self.ga_config["sizeFac"]}, 
                    popSize={self.ga_config["popSize"]},  elitism={self.ga_config["elitism"]}, 
                    stallGenMax={self.ga_config["stallGenMax"]}, relTol={self.ga_config["relTol"]},
                    pMutation={self.ga_config["pMutation"]}, selPress={self.ga_config["selPress"]},
                    NAFac={self.ga_config["NAFac"]},
                    maxTime={self.ga_config["maxTime"]}, verbose={self.ga_config["verbose"]})
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
                t <- system.time(resILP <- CellNOptR:::createAndRunILP(
                    model, md, cnolistReal, accountForModelSize = {self.ilp_config['accountForModelSize']}, 
                    sizeFac = {self.ilp_config['sizeFac']}, mipGap={self.ilp_config['mipGap']}, 
                    relGap={self.ilp_config['relGap']}, timelimit={self.ilp_config['timelimit']}, 
                    cplexPath = "{self.ilp_config['cplexPath']}", method = "{self.ilp_config['method']}", 
                    numSolutions = {self.ilp_config['numSolutions']}, limitPop = {self.ilp_config['limitPop']}, 
                    poolIntensity = {self.ilp_config['poolIntensity']}, poolReplace = {self.ilp_config['poolReplace']})
                    );
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
            optModel <- processOptimizedModel(optModel)
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
            # SIFToBoolNet(sifFile     = "{sif_fname}",
            #             boolnetFile = "{boolnet_fname}",
            #             CNOlist     = cnolist,
            #             model       = optModel,
            #             fixInputs   = FALSE,
            #             preprocess  = TRUE,
            #             ignoreAnds  = TRUE)
            
            result <- writeBnetFromModel(optModel, "{boolnet_fname}")

            verifyBoolNetConversion(optModel, "{boolnet_fname}")
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
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(output_dir, fname)
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")
        
        AA = AttractorAnalysis(self.GD_MODEL, opt_fname)
        results = AA.comparison()
        
        if 'ilp' in output_dir:
            results['method']     = "ILP"
        else:
            results['method']     = "GA"
        # print("Total time taken for evaluation:", r('t[["elapsed"]]'))
        results['total_time']         = limit_float(r('t[["elapsed"]]'))
        results['change_percent']     = limit_float(self.ChangePct)
        # results.to_csv(os.path.join(output_dir, "results.csv"), index=False)

        return results  
    
    def run_full_analysis(self, method="ga"):
        # Load network and data
        model = self.load_network()
        cnolist = self.load_data()
        
        # Preprocess
        self.preprocess_network()

        # Optimize
        output_file = os.path.join(self.output_file, method)                
        os.makedirs(output_file, exist_ok=True)
        print(f"Optimizing with method: {method}")
        opt_results = self.optimize_network(PLOT=True, method=method, output_file=output_file)

        # Save results
        self.save_results(output_file)
        print(f"Results saved to {output_file}/")
        
        # Evaluate and cross-validate
        print("Evaluating and cross-validating the model...")
        results = self.evaluate_model(output_file)        
        print("Analysis completed successfully!")
        return model, cnolist, results

def main(dataset="toy", method="ga", manager=None):
    """Example usage of CellNOptAnalyzer"""
    
    # file = "data/apoptosis.xml"
    # For MIDAS file, you'll need to create or provide one
    # This is just an example - you'll need actual experimental data
    # midas_file = "data/experimental_data.csv"  # You need to create this

    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer(
        dataset=dataset,
        manager=manager
    )

    model, cnolist, results = analyzer.run_full_analysis(
        method=method
    )
    return model, cnolist, results

if __name__ == "__main__":
    # model, cnolist, results = main(Test=True, method="ilp")
    from tools.config import NetworkPerturbationConfig, AdaptiveParameterManager
    config = NetworkPerturbationConfig(
        change_percent=0.,
        size_adaptation_strength=2.0,
        generalization_focus=True
    )
    manager = AdaptiveParameterManager(config)
    model, cnolist, results = main(dataset="toy", method="ga", manager=manager)