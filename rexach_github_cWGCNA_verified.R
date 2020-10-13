
## This script takes four expression datasetes and creates 1. consensus networks and 2. consensus gene module connectivity scores (kME)

R-3.0.1

rm(list=ls())

#require('WGCNA')

library(WGCNA) # version  1.27.1 

enableWGCNAThreads()

options(StringsAsFactors=F)



## set working directory to code directory

# load normalized TPR50 cortex RNAseq data (outliers removed) (published - Swarup et al, 2019) and cleaned and normalized Tg4510 microglia RNAseq data (published - Wang et al., 2018); together with sample metadata


load("input_datExpr.rda")

load("input_targets.rda")


## subset dataset to common genes 


gnS=intersect(colnames(datExpr_micro),intersect(colnames(datExpr_C57BL6),intersect(colnames(datExpr_DBA),colnames(datExpr_FVB))))
length(gnS)

datExpr_micro= datExpr_micro[,na.omit(match(gnS,colnames(datExpr_micro)))]

datExpr_C57BL6= datExpr_C57BL6[,na.omit(match(gnS,colnames(datExpr_C57BL6)))]

datExpr_FVB= datExpr_FVB[,na.omit(match(gnS,colnames(datExpr_FVB)))]

datExpr_DBA= datExpr_DBA[,na.omit(match(gnS,colnames(datExpr_DBA)))]


# format combined data for consensus WGCNA

nSets=4

setLabels=c("Tg4510_microglia","TPR50.DBA","TPR50.FVB","TPR50.C57")

shortLabels=c("Tg4510_microglia","TPR50.DBA","TPR50.FVB","TPR50.C57")

multiExpr=vector(mode="list",length=nSets)

multiExpr[[1]] = list(data=as.data.frame(datExpr_micro)) 

names(multiExpr[[1]]$data)=colnames(datExpr_micro)

rownames(multiExpr[[1]]$data)=rownames(datExpr_micro)


multiExpr[[2]] = list(data=as.data.frame(datExpr_DBA)) 

names(multiExpr[[2]]$data)=colnames(datExpr_DBA)

rownames(multiExpr[[2]]$data)=rownames(datExpr_DBA)


multiExpr[[3]] = list(data=as.data.frame(datExpr_FVB)) 

names(multiExpr[[3]]$data)=colnames(datExpr_FVB)

rownames(multiExpr[[3]]$data)=rownames(datExpr_FVB)


multiExpr[[4]] = list(data=as.data.frame(datExpr_C57BL6)) 

names(multiExpr[[4]]$data)=colnames(datExpr_C57BL6)

rownames(multiExpr[[4]]$data)=rownames(datExpr_C57BL6)


checkSets(multiExpr) # check data size and content



# format metadata for consensus WGCNA 



multiMeta=list(Tg4510_microglia

 =list(data=targets_micro),TPR50.DBA=list(data=targets_DBA),TPR50.FVB=list(data=targets_FVB),TPR50.C57=list(data=targets_C57BL6))

checkSets(multiMeta) # check data content


## network construction


# test and chose soft-thresholding powers to achieve scale-free topology

powers = c(seq(1,10,by=1), seq(12,40, by=2));

# initialize a list to hold the results of scale-free analysis

powerTables = vector(mode = "list", length = nSets);

# call the network topology analysis function for each set in turn

for (set in 1:nSets)

	powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
	
	verbose = 5,networkType="signed",corFnc="bicor")[[2]]);


# plot the results:

pdf("../figures/softPower.pdf", height=10, width=18)

colors = c("blue", "red","black")

plotCols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",

"Max connectivity");


# get the minima and maxima of the plotted points

ylim = matrix(NA, nrow = 2, ncol = 4);

	for (set in 1:nSets)

	{

	for (col in 1:length(plotCols))

		{

			ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);

			ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);

		}

	}


# plot the quantities in the chosen columns vs. the soft thresholding power

par(mfcol = c(2,2));

par(mar = c(4.2, 4.2 , 2.2, 0.5))

cex1 = 0.7;

	for (col in 1:length(plotCols)) for (set in 1:nSets)

	{

		if (set==1)

		{

			plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],

			xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],

			main = colNames[col]);

			addGrid();

		}

		if (col==1)

		{

		text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],

		labels=powers,cex=cex1,col=colors[set]);

		} else

		text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],

		labels=powers,cex=cex1,col=colors[set]);

		if (col==1)

		{

		legend("bottomright", legend = setLabels, col = colors, pch = 20) ;

		} else

		legend("topright", legend = setLabels, col = colors, pch = 20) ;

	}

dev.off()





## construct consensus WGCNA network 


softPower=14


net=blockwiseConsensusModules(multiExpr, blocks = NULL, 

                                         maxBlockSize = 30000,
                                         
                                         randomSeed = 12345,
                                         
                                         corType = "pearson",
                                         
                                         consensusQuantile=0, 
                                         
                                         power = softPower,
                                         
                                         networkType = "signed",
                                         
                                         TOMType = "unsigned",
                                         
                                         TOMDenom = "min",
                                         
                                         scaleTOMs = TRUE, scaleQuantile = 0.8,
                                         
                                         sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                         
                                         useDiskCache = TRUE, chunkSize = NULL,
                                         
                                         deepSplit = 2,
                                         
                                         detectCutHeight = 0.995, minModuleSize = 100,
                                         
                                         mergeCutHeight = 0.2,
                                         
                                         saveTOMs = TRUE,
                                         
                                         consensusTOMFileNames = "ConsensusTOM-block.%b.rda")
                                         
                                         

consMEs = net$multiMEs;

moduleLabels = net$colors;


# convert the numeric labels to color labels

moduleColors = labels2colors(moduleLabels)

consTree = net$dendrograms[[1]];

pdf("SignedDendro_Consensus2.pdf",height=10, width=15)

plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,

main = "Consensus gene dendrogram and module colors")

dev.off()
 


load("ConsensusTOM-block.1.rda")


# various tree cutting parameter

    consTree= hclust(1-consTomDS,method="average");






mColorh <- mLabelh <- colorLabels <- NULL

  for (minModSize in c(40,100)) {
  	
    for (dthresh in c(0.2,0.25)) {
    	
      for (ds in c(2,4)) {
      	
          print("Trying parameters:")
          
          print(c(minModSize,dthresh,ds))
          
          tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
          
            minClusterSize = minModSize, cutHeight = 0.995,
            
            deepSplit = ds, distM = as.matrix(1-consTomDS))
            
          
          merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels, cutHeight = dthresh)
                                      
          mColorh <- cbind(mColorh,labels2colors(merged$colors))
          
          mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
          
        }
      }
    }



### first network - microglia

tmpMulti =vector(mode="list",length=nSets)

thisMeta <- multiMeta[[1]]$data

thisExpr <- multiExpr[[1]]$data
		
datTraits<- targets_micro[,c(2:4)]
  
cond1=as.character(datTraits$condition)

ind=which(cond1=="Tg")

cond1[ind]="Tg"

condition1 =as.numeric(factor(cond1, c("WT","Tg"))) -1 

		
		  tmpMulti[[1]]$traitmat <- as.data.frame(cbind(as.character(datTraits[,1]),as.numeric(datTraits[,2]),condition1))
		  
		  rownames(tmpMulti[[1]]$traitmat) <- rownames(datTraits)
		  
		  colnames(tmpMulti[[1]]$traitmat) <- c("MouseID","age","condition")
		  
		
		  geneSigs <- matrix(NA,nrow=1,ncol=ncol(thisExpr))
		  
		  ## Find adjusted multiple R^2 for each gene withe ach categorical variable
		  
		  for (i in 1:ncol(geneSigs)) {
		  	
		    exprvec <- as.numeric(thisExpr[,i])
		    
			conditionr=cor(tmpMulti[[1]]$traitmat[,3],exprvec) 
	
			geneSigs[,i]=c(conditionr)
		  }
		  
		  tmpMulti[[1]]$genesigs <- geneSigs
		
		  geneSigs <- numbers2colors(as.numeric(geneSigs),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
		  
	      colnames(geneSigs) <- c("TgCondition.microglia")
	      
		  tmpMulti[[1]]$genecols <- t(geneSigs)
		  
		  tmpMulti[[1]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0 and power of 24"))
		  
		  tmpMulti[[1]]$netData$TOMdendrogram <- consTree
		  
		  tmpMulti[[1]]$netData$moduleColors <- mColorh
		  
		  tmpMulti[[1]]$netData$cutParameters <- mLabelh
		  
		  tmpMulti[[1]]$netData$annotColors <- geneSigs
		  


### second network -  DBA
  
thisMeta <- multiMeta[[2]]$data
 
thisExpr <- t(multiExpr[[2]]$data)

tmpMulti[[2]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"Sample.ID"]))

rownames(tmpMulti[[2]]$traitmat) <- rownames(thisMeta)
  
colnames(tmpMulti[[2]]$traitmat) <- c("Wt.Tg","Sample.ID")

geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  
 
 ## find adjusted multiple R^2 for each gene with each categorical variable
  
 	 for (i in 1:ncol(geneSigs)) {
  	
   		 exprvec <- as.numeric(thisExpr[i,])
    
   		conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[2]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    
    		rinr <- bicor(tmpMulti[[2]]$traitmat[,"Sample.ID"],exprvec)
    
    		geneSigs[,i] <- c(conditionr,rinr)
    
 	 	}
  
		tmpMulti[[2]]$genesigs <- geneSigs
  
 		geneSigs[1,] =  numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
   
		geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
   
		rownames(geneSigs) <- c("DBA_Wt.Tg","DBA_sampleid")
  
  		tmpMulti[[2]]$genecols <- geneSigs

  		tmpMulti[[2]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  
  		tmpMulti[[2]]$netData$TOMdendrogram <- geneTree
  
  		tmpMulti[[2]]$netData$moduleColors <- mColorh
  
  		tmpMulti[[2]]$netData$cutParameters <- mLabelh
  
 		tmpMulti[[2]]$netData$annotColors <- geneSigs


### third network - FVB

thisMeta <- multiMeta[[3]]$data
 
thisExpr <- t(multiExpr[[3]]$data)

tmpMulti[[3]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"Sample.ID"]))

rownames(tmpMulti[[3]]$traitmat) <- rownames(thisMeta)
  
colnames(tmpMulti[[3]]$traitmat) <- c("Wt.Tg","Sample.ID")

geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  
  ## Find adjusted multiple R^2 for each gene with each categorical variable
  
  	for (i in 1:ncol(geneSigs)) {
  	
    		exprvec <- as.numeric(thisExpr[i,])
    
    		conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[3]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    
    		rinr <- bicor(tmpMulti[[3]]$traitmat[,"Sample.ID"],exprvec)
    
    		geneSigs[,i] <- c(conditionr,rinr)
    
  		}
  
  		tmpMulti[[3]]$genesigs <- geneSigs

  		geneSigs[1,] =  numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  
  		geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  
  		rownames(geneSigs) <- c("FVB_Wt.Tg","FVB_Sample.ID")
  
  		tmpMulti[[3]]$genecols <- geneSigs

  		tmpMulti[[3]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  
  		tmpMulti[[3]]$netData$TOMdendrogram <- geneTree
  
  		tmpMulti[[3]]$netData$moduleColors <- mColorh
  
  		tmpMulti[[3]]$netData$cutParameters <- mLabelh
  
  		tmpMulti[[3]]$netData$annotColors <- geneSigs


### fourth network - C57

thisMeta <- multiMeta[[4]]$data

thisExpr <- t(multiExpr[[4]]$data)

tmpMulti[[4]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"Sample.ID"]))

rownames(tmpMulti[[4]]$traitmat) <- rownames(thisMeta)
  
colnames(tmpMulti[[4]]$traitmat) <- c("Wt.Tg","Sample.ID")

geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  
  ## Find adjusted multiple R^2 for each gene with each categorical variable
  
	for (i in 1:ncol(geneSigs)) {
  	
    		exprvec <- as.numeric(thisExpr[i,])
    
    		conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[4]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    
    		rinr <- bicor(tmpMulti[[4]]$traitmat[,"Sample.ID"],exprvec)
    
    		geneSigs[,i] <- c(conditionr,rinr)
    
  		}
  
  		tmpMulti[[4]]$genesigs <- geneSigs
  
 		geneSigs =   tmpMulti[[4]]$genesigs

  		geneSigs[1,] <- numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  
  		geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  
  		rownames(geneSigs) <- c("C57_Wt.Tg","C57_sampleid")
  
  		tmpMulti[[4]]$genecols <- geneSigs

  		tmpMulti[[4]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  
  		tmpMulti[[4]]$netData$TOMdendrogram <- geneTree
  
  		tmpMulti[[4]]$netData$moduleColors <- mColorh
  
  		tmpMulti[[4]]$netData$cutParameters <- mLabelh
  
  		tmpMulti[[4]]$netData$annotColors <- geneSigs
  
  
###### combine data	and plot final dendrogram
		
	
		  mColorh1 <- cbind(mColorh,t(tmpMulti[[1]]$genecols),t(tmpMulti[[2]]$genecols),t(tmpMulti[[3]]$genecols),t(tmpMulti[[4]]
		  $genecols))
		  
		  mLabelh1 <- c(mLabelh,rownames(tmpMulti[[1]]$genecols),rownames(tmpMulti[[2]]$genecols),rownames(tmpMulti[[3]]
		  $genecols),rownames(tmpMulti[[4]]$genecols))
		  
		  pdf("Dendrogram_with_GeneSigs.pdf",height=30,width=25)
		  
		  plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Microglia_TPR50 Data")
		  
		  multiData.Tg4510_microglia <- tmpMulti[[1]]
		  
		  multiData.TPR50_DBA <- tmpMulti[[2]]
		  
		  multiData.TPR50_FVB <- tmpMulti[[3]]
		  
		  multiData.TPR50_C57 <- tmpMulti[[4]]

		  dev.off()




mms=100

ds=4

dthresh=0.25


tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,

       	minClusterSize = mms, cutHeight = 0.995, deepSplit = ds, distM = as.matrix(1-consTomDS))
              
       	merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,cutHeight = dthresh)     
         
       	mColorh <- cbind(labels2colors(merged$colors),t(multiData.Tg4510_microglia$genecols), 
        
       	t(multiData.TPR50_DBA$genecols),t(multiData.TPR50_FVB$genecols),t(multiData.TPR50_C57$genecols))
        
       	mLabelh <- c("Merged 
       		Colors",rownames(multiData.Tg4510_microglia$genecols),rownames(multiData.TPR50_FVB$genecols),rownames(multiData.TPR50_DBA$genecols),rowname

		s(multiData.TPR50_C57$genecols))

 
 
  pdf("ConsensusTOM_withGeneSigs.pdf",height=15,width=25)
  
  plotDendroAndColors(consTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with 
  
  power =30, mms=",mms,"ds=",ds,"dthresh=",dthresh));
  
  dev.off()







#### secondary clustering of genes from selected modules

# this code is to test for additional associations among genes assigned to two different but biologically-related modules 

# load input data

datExpr =rbind(datExpr_micro, datExpr_C57BL6, datExpr_FVB, datExpr_DBA) # combined expression data

geneInfo = read.csv("consensusKMEtable") # consensus network kME table



# extract modules of interest for secondary clustering

turq = geneInfo[which(geneInfo[,3]=="turquoise"),] 

gy = geneInfo[which(geneInfo[,3]=="greenyellow"),]

gyturq = rbind(turq,gy)

modgenes=gyturq$Ensembl.Gene.ID

geneNames=gyturq$GeneSymbol

gnS=intersect(modgenes,colnames(datExpr))

thisExpr=datExpr[,match(gnS,colnames(datExpr))]
 
gyt=gyturq[match(gnS,modgenes),]

colnames(thisExpr) = toupper(gyt$GeneSymbol)

# subset for genes with experimentally annotated PPI (this option was selected to hone in on genes in signaling pathways)

load("./BGandIWcombinedPPI_5-8-2014.Rdata")

keepgenes = intersect(rownames(ppiMat),colnames(thisExpr))

coexpInd = na.omit(match(keepgenes,colnames(thisExpr)))

ppiInd = na.omit(match(keepgenes,rownames(ppiMat)))

seedExpr = thisExpr[,coexpInd]


# generate bicor matrix 

adjMat = bicor(seedExpr) 


# remove genes with negative bicor prior to TOM

adjMat[adjMat < 0] <- NA 


# calculate topological overlap
  
TOM = TOMsimilarity(adjMat)

dissTOM = 1-TOM

geneTree = flashClust(as.dist(dissTOM), method = "average");


# assign genes to new secondary modules


mms=40

ds =2

dthresh=0.1
		
library(flashClust)
		
geneTree = flashClust(as.dist(dissTOM), method = "average");

tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
        
merge <- mergeCloseModules(exprData = seedExpr,colors = tree$labels, cutHeight = dthresh)

mergedColors = labels2colors(merge$colors);

mergedMEs = merge$newMEs;

MEList=moduleEigengenes(seedExpr,colors = mergedColors,softPower= 7, nPC=1)

MEs=MEList$eigengenes

MEs=orderMEs(MEs)

moduleColors = mergedColors

mColorh <- cbind(labels2colors(merge$colors))

mLabelh <- c("Merged Colors")
        

pdf("secondary_modules.pdf",height=10,width=16) 

	plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGui0de=TRUE,dendroLabels=FALSE,main= paste("Signed cor network with 
	
	power 4","mms=",mms,"ds=",ds,"dthresh=",dthresh));

	dev.off()
  
  
  
# calculate gene-module connectivity scores


KMEs<-signedKME(seedExpr, MEs,outputColumnName = "kME",corFnc = "bicor")


          













