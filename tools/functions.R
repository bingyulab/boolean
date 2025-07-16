library(CellNOptR)


write_boolnet_to_sif <- function(net, file) {
  sif <- data.frame()
  for (gene in names(net$interactions)) {
    expr <- net$interactions[[gene]]$expression
    # Remove parentheses and split by operators
    expr_clean <- gsub("[()]", "", expr)
    # Split by | or & (OR/AND), then trim whitespace
    regulators <- unlist(strsplit(expr_clean, "[|&]"))
    regulators <- trimws(regulators)
    regulators <- regulators[regulators != "" & regulators != gene]
    for (reg in regulators) {
      if (startsWith(reg, "!")) {
        interaction <- -1
        reg_clean <- substring(reg, 2)
      } else {
        interaction <- 1
        reg_clean <- reg
      }
      sif <- rbind(sif, data.frame(source=reg_clean, interaction=interaction, target=gene, stringsAsFactors=FALSE))
    }
  }
  write.table(sif, file=file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Usage:
# write_boolnet_to_sif(net, "network.sif")

SIFToBoolNet <- function(sifFile, boolnetFile, CNOlist, model=NULL, fixInputs=TRUE, preprocess=TRUE, ignoreAnds=TRUE)
{
  if (preprocess)
  {
    if (is.null(model))
      sif <- readSIF(sifFile)
    else
      sif = model
    sif <- preprocessing(data=CNOlist, model=sif)
    system("rm tmp.sif")
    writeSIF(sif,file="tmp.sif", overwrite=TRUE)
    sifFile <- "tmp.sif"
  }
    
  if ((class(CNOlist)=="CNOlist")==FALSE){
      CNOlist = CellNOptR::CNOlist(CNOlist)
  } 

  sif <- read.table(sifFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
  if (ncol(sif) != 3)
    sif <- read.table(sifFile,sep=" ",header=FALSE,stringsAsFactors=FALSE)
  
  genes <- unique(c(sif[,1],sif[,3],colnames(CNOlist@cues)))
  genes <- gsub("-","_",genes,fixed=TRUE)
  sif[,1] <- gsub("-","_",sif[,1],fixed=TRUE)
  sif[,3] <- gsub("-","_",sif[,3],fixed=TRUE)  
  #missingGenes <- sort(setdiff(genes,unique(c(cues,signals))))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following network genes were not present in the MIDAS file:",paste(missingGenes,collapse=",")))
  #missingGenes <- sort(setdiff(unique(c(cues,signals)),genes))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following MIDAS genes were not present in the network:",paste(missingGenes,collapse=",")))
  
  #cat("Both cues and signals:\n")
  #print(intersect(cues,signals))
  
  andIdx <- grep("and",genes,fixed=TRUE)
  andIdx <- andIdx[order(genes[andIdx],decreasing=TRUE)]
  ands <- genes[andIdx]

  andStrings <- sapply(ands, function(gene)
  {
    facIdx <- which(sif[,3] == gene)
    if (ignoreAnds && length(facIdx) > 1) 
      return(NA)
      
    if (length(facIdx) > 0)
      paste("(",paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" & "),")",sep="")
  })
  
  if (ignoreAnds)
    ignoredAnds <- ands[is.na(andStrings)]
  else
    ignoredAnds <- c()
    
  geneStrings <- sapply(genes[-andIdx], function(gene)
  {
    facIdx <- which(sif[,3] == gene & !(sif[,1] %in% ignoredAnds))    
      
    if (length(facIdx) > 0)
      paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" | ")
    else
    if (fixInputs && (gene %in% colnames(CNOlist@cues)) && 
        !(gene %in% colnames(CNOlist@signals[[1]])))
      paste(gene, "[1.0]", sep="")
    else    
      gene
  })
  
  geneStrings <- sapply(geneStrings, function(s)
  {
    
    for (a in 1:length(andIdx))
    {
      s <- gsub(ands[a], andStrings[a], s, fixed=TRUE)
    }
    s
  })
  genes <- genes[-andIdx]
  print(genes)
  #print(outputStrings)
  sink(boolnetFile)
  cat("targets, factors\n")  
  for (i in 1:length(genes))
  {
    cat(genes[i],", ",geneStrings[i],"\n",sep="")
  }    
  sink()
}


# library(CellNOptR)
# CNOl <- makeCNOlist(readMIDAS(midasFile),subfield=FALSE)
# if ((class(CNOl)=="CNOlist")==FALSE){
#     CNOl = CellNOptR::CNOlist(CNOl)
# }
# CNOl <- makeCNOlist(readMIDAS(midasFile),subfield=FALSE)
# SIFToBoolNet(sifFile, "model_draft.txt", CNOl, fixInputs=fixInputs, preprocess=preprocess)
  

simulate_CNO <- function(model, bString, simList=NULL, CNOlist, indexList=NULL)
{

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
    if (is.null(simList)==TRUE){
        simList <- prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # keep simList and indxList for back compatibility ?
    modelCut <- cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim0[,indexList$inhibited] <- 1 - Sim0[,indexList$inhibited]
    simRes0 <- as.matrix(Sim0[,indexList$signals])

    # t1
    Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim[,indexList$inhibited] <- 1 - Sim[,indexList$inhibited]
    
    colnames(Sim) <- model$namesSpecies
    simRes <- as.matrix(Sim[,indexList$signals])

    sig <- #as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
              as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
    
    simResults <- list(input=CNOlist@cues,
                       t0=simRes0, t1=simRes, trueSig=CNOlist@signals[[2]])

    return(simResults)
}


pkn2sif<-function(model,optimRes=NA,writeSif=FALSE, filename="Model"){
	
  if (is.na(optimRes[1])){
    BStimes<-rep(1,length(model$reacID))
  }else{
    BStimes<-optimRes$bString
  }	

	findOutput<-function(x){
		sp<-which(x == 1)
		sp<-model$namesSpecies[sp]
		}
		
	reacOutput<-apply(model$interMat,2,findOutput)
	
	findInput<-function(x){
		sp<-which(x == -1)
		sp<-model$namesSpecies[sp]
		}
		
	reacInput<-apply(model$interMat,2,findInput)
		
#if the class of reacInput is not a list, then there are no AND gates
	if(!is(reacInput,"list")){
	
		isNeg<-function(x){
			isNegI<-any(x == 1)
			return(isNegI)
			}
		
    if (all(model$notMat == 0)) {
      # No negative connections: all are positive
      inpSign <- rep(1, ncol(model$interMat))
    } else {
      inpSign <- apply(model$notMat, 2, isNeg)
      inpSign <- !inpSign
      inpSign[inpSign] <- 1
      inpSign[!inpSign] <- -1
    }
    
		sifFile<-cbind(reacInput,inpSign,reacOutput)
    sifFile<-sifFile[BStimes==1,]
		colnames(sifFile)=NULL
    rownames(sifFile)=NULL
    
		}else{
		
#in this case there are AND gates and so we need to create dummy "and#" nodes
			sifFile<-matrix(ncol=3)
			nANDs<-1
			for(i in 1:length(reacOutput)){
			  if (BStimes[i]==1){

				  if(length(reacInput[[i]]) == 1){
            tmp<-matrix(0,nrow=1,ncol=3)
					  tmp[1,1]<-reacInput[[i]]
					  tmp[1,3]<-reacOutput[i]
					  tmp[1,2]<-ifelse(
						  any(model$notMat[,i] == 1),-1,1)
              sifFile<-rbind(sifFile,tmp)
          
					}else{
					
						for(inp in 1:length(reacInput[[i]])){
              tmp<-matrix(0,nrow=1,ncol=3)
							tmp[1,1]<-reacInput[[i]][inp]
							tmp[1,3]<-paste("and",nANDs,sep="")
							tmp[1,2]<-ifelse(
								model$notMat[which(reacInput[[i]][inp]==rownames(model$notMat)),i] == 1,
								-1,1)
              sifFile<-rbind(sifFile,tmp)
							}
						tmp<-matrix(0,nrow=1,ncol=3)	
						tmp[1,1]<-paste("and",nANDs,sep="")
						tmp[1,3]<-reacOutput[i]
						tmp[1,2]<-1
            sifFile<-rbind(sifFile,tmp)
						
            nANDs<-nANDs+1
            }
			    }
				}
				
      sifFile<-sifFile[2:dim(sifFile)[1],]
			
			}

#this is the edge attribute matrix that contains, for each edge, whether it is
#absent from the model (0), present at t1(1) or present at t2(2)
  if (writeSif==TRUE){
	filename<-paste(filename, ".sif", sep="")
    writeSIF(sifFile, filename=filename)
  }

  return(sifFile)
}

toSIF <- function(model, filename, overwrite=FALSE){

    # convert internal model structure to a SIF matrix
    sif = pkn2sif(model)

    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(sif[,1:3],file=filename,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
    else{
       if (overwrite==TRUE){
            write.table(sif[,1:3],file=filename,
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
        }
        else{
           stop(paste("File ", filename, "already exists.",  sep=""))
        }
    }


}