library(CellNOptR)


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
    writeSIF(sif,file="tmp.sif")
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
  