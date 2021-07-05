#### SCENIC - module of expression ####
SCEN<-function(){
	options(max.print=50)
	table(GEX.backup$Condition, GEX.backup$CellTypeFinal)
	# NOTE: only 19 CD4 T cell and 0 CD8 T cell in Naive condition: not used
	MAIN_WD=getwd()
	
	# The correlation has to be run per condition, otherwise we mix everything
	### Load data
	for( index in 2){
		if(index==1){ setwd(MAIN_WD) ; GEX<- subset(GEX.backup, Condition=="Peak") ; SUBDIR="SCENIC.Peak" ; dir.create(SUBDIR) ; setwd(SUBDIR) }
		if(index==2){ setwd(MAIN_WD) ; GEX <- subset(GEX.backup, Condition=="Remission") ; SUBDIR="SCENIC.Recovery"; dir.create(SUBDIR); setwd(SUBDIR) }

		exprMat <- as.matrix(GetAssayData(GEX, "counts")) #scale.data counts
		cellInfo <- as.data.frame(GEX@meta.data)
		cellInfo$nGene <- colSums(exprMat>0)
		rownames(exprMat) # this should be genes #exprMat <- t(exprMat)
	
		### Initialize settings
		scenicOptions <- initializeScenic(org="mgi",dbDir="../SCENIC.feather",datasetTitle=SUBDIR,nCores=4) # we can ignore the missing column errorhttps://github.com/aertslab/SCENIC/issues/168
		saveRDS(scenicOptions, file="int/scenicOptions.Rds")
	
		### Co-expression network
		genesKept <- SCENIC::geneFiltering(exprMat, scenicOptions) # save in 1.1_genesKept.Rds   #it seem that this has to use the raw counts, since we convert to log2 later
		exprMat_filtered <- exprMat[genesKept, ]
		SCENIC::runCorrelation(exprMat_filtered, scenicOptions) # save in 1.2_corrMat.Rds
		exprMat_filtered_log <- log2(exprMat_filtered+1) 
		SCENIC::runGenie3(exprMat_filtered_log, scenicOptions)
		
		### Build and score the GRN - gene regulatory network
		exprMat_log <- log2(exprMat+1)
		#Step 1: Convert the output from GENIE3/GRNBoost to co-expression modules 
		scenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
		#Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis) 
		scenicOptions <- SCENIC::runSCENIC_2_createRegulons(scenicOptions, minGenes = 10) # , coexMethod=c("top5perTarget")
		#Step 3: AUCell (scoring the regulons on the individual cells) 
		scenicOptions <- SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMat_log) # why exprMat_log in tuto ??
		
		saveRDS(scenicOptions, file="int/scenicOptions.Rds")
		setwd(MAIN_WD)
		rm(scenicOptions,exprMat_log,exprMat_filtered_log,exprMat_filtered,genesKept,example_avgexp_data,exprMat,cellInfo)
	}
	
}
