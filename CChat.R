#### CellChat for CCC - cell cell communication / CCI cell cell interactions ####
# Jin et al. (2020)
CChat<-function(){
	options(stringsAsFactors = FALSE) # !!!!
	NAMEIDENT="CellTypeFinal" #"annot_RNA_snn_res.0.1_CIPR" 
	for( ORIGIN in unique( GEX$Condition ) ) {
		
		GEX <- subset( GEX.backup, Condition == ORIGIN)
		GEX$CellTypeFinal<-factor(GEX$CellTypeFinal)
		
		for( SELECT in c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor")){
			OUTDIR=paste0('15_cellchat_',ORIGIN,"_",gsub(" ","",SELECT)) ; dir.create(OUTDIR,showWarnings = F)
			print(paste(">>>>>>>>> Currently:",OUTDIR))
			### setup
			DATABASE="CellChatDB.mouse"
			PPI=PPI.mouse
			
			# setup object
			# in the future, when no bugs
			#cellchat <- createCellChat(object = GEX)
			
			Idents(GEX)<-NAMEIDENT
			cellchat <- createCellChat(object = GetAssayData(GEX,assay="RNA",slot="data"),meta=GEX@meta.data,group.by=NAMEIDENT)
			cellchat@DB <- subsetDB(CellChatDB.mouse,search=SELECT) # use Secreted Signaling for cell-cell communication analysis
			
			# process (warning: use the proper input data ~ species: PPI.mouse)
			if(T){
				### We first identify over-expressed ligands or receptors in one cell group, 
				#and then project gene expression data onto protein-protein interaction (PPI) network.
				#The over-expressed ligand-receptor interactions are identified if either the ligand or receptor is over-expressed.
				
				cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
				#future::plan("multiprocess", workers = 4) # do parallel
				cellchat <- identifyOverExpressedGenes(cellchat)	
				cellchat <- identifyOverExpressedInteractions(cellchat)
				cellchat <- projectData(cellchat,PPI)
				
				
				#Inference of cell-cell communication network
				cellchat <- computeCommunProb(cellchat) #,raw.use=F,trim=0.1,type="truncatedMean",population.size=T)
				
				# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
				cellchat <- filterCommunication(cellchat,min.cells=10)
				
				#Infer the cell-cell communication at a signaling pathway level
				cellchat <- computeCommunProbPathway(cellchat)
				
				#Calculate the aggregated cell-cell communication network
				cellchat <- aggregateNet(cellchat)
				
				saveRDS(cellchat, paste0(OUTDIR,"/cellchat.mouse.Rds"))
			}
			#extract signif interaction - save in excel
			if(T){
				#cellchat@meta
				#returns a data frame consisting of all the inferred cell-cell communications at the level of 
				#ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the 
				#level of signaling pathways
				rm(df.net)
				try({ df.net <- subsetCommunication(cellchat) })
				if( exists("df.net") ) { write.xlsx(df.net,paste0(OUTDIR,"/","df.net.xlsx")) }else{print("no df.net")}
				
			}
		
			rm(df.net,cellchat)
			try({ df.net<-read.xlsx(paste0(OUTDIR,"/","df.net.xlsx")) })
			
			if ( exists("df.net") ){ #& nrow(df.net) > 0
				### If already processed:
				cellchat<-readRDS(paste0(OUTDIR,"/cellchat.mouse.Rds"))
				
				#Visualization and systems analysis of cell-cell communication network
				
				# circular links: not amazing since all at once
				if(T){
					OUTNAME=paste0(OUTDIR,"/",DATABASE,".netVisual_circle.number.png");png(OUTNAME,width=800,height=800)
					netVisual_circle(cellchat@net$count,vertex.weight=as.numeric(table(cellchat@idents)),weight.scale=T,label.edge=F,title.name="Number of interactions")
					dev.off()
					
					OUTNAME=paste0(OUTDIR,"/",DATABASE,".netVisual_circle.strength.png");png(OUTNAME,width=800,height=800)
					netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale=T,label.edge=F,title.name="Interaction weights/strength")
					dev.off()
					
					rm(OUTNAME)
				}
				
				# circular links SPLIT per celltype // circle // # This split the previous "weighted" circular plot
				#summarizing the communication probability
				if(T){
					dir.create(paste0(OUTDIR,"/circle/"))
					groupSize <- as.numeric(table(cellchat@idents))
					mat <- cellchat@net$weight
					for (index in 1:nrow(mat)) {
						#create new emtpy matrix
						mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
						#fill the new patrix with only data from the current index/ celltype
						mat2[index, ] <- mat[index, ]
						NAME=rownames(mat)[index]
						NAME=gsub("[/]","-",NAME)
						print(paste("index",index,":",NAME))
						OUTNAME=paste0(OUTDIR,"/circle/",NAME,".pdf");pdf(OUTNAME,width=8,height=8,useDingbats=F)
						netVisual_circle(mat2,vertex.weight=groupSize,
										 weight.scale=T,
										 edge.weight.max=max(mat),
										 title.name=rownames(mat)[index],vertex.label.cex = .8)
						dev.off()
					}
					rm(mat,mat2,NAME,OUTNAME,mat,groupSize,index)
				}
				
				# Compute and visualize the contribution of each (significant)
				# ligand-receptor pair in the overall signaling pathways 
				if(T){
					SIGNALING=cellchat@netP$pathways
					
					WIDTH=10
					HEIGHT=10
					FONT_SIZE=12
					
					OUTNAME=paste0(OUTDIR,"/","L-R pair Contribution.pdf")
					pdf(OUTNAME,width=WIDTH,height=HEIGHT)
					CONTRIB=netAnalysis_contribution(cellchat,signaling=SIGNALING,thresh=0.05,return.data=T,font.size=FONT_SIZE)
					print(CONTRIB$gg.obj)
					dev.off()
					
					#CONTRIB$LR.contribution$name
					CONTRIB$LR.contribution$name = factor(CONTRIB$LR.contribution$name,levels=CONTRIB$LR.contribution$name[order(CONTRIB$LR.contribution$contribution)])
					CONTRIB$LR.contribution <- subset( CONTRIB$LR.contribution[ CONTRIB$LR.contribution$name != "1",])
					PLOT=ggplot(CONTRIB$LR.contribution) +
						geom_bar(aes(x=contribution*100,y=name),stat="identity") +
						theme_minimal(base_size=FONT_SIZE) +
						ggtitle("Contribution of each L-R pair\nin the overall signaling pathways") +
						xlab("Relative contribution (%)")+
						ylab("") 
					
					OUTNAME=paste0(OUTDIR,"/","L-R pair Contribution.2.pdf");pdf(OUTNAME,width=WIDTH,height=HEIGHT);print(PLOT);dev.off()
					#sum(CONTRIB$LR.contribution$contribution)
					rm(CONTRIB,PLOT,OUTNAME)
				}
				
				# Bubble plot
				if(T){
					dir.create( paste0(OUTDIR,"/","BubblePlot/" ) )
					for( SOURCENAME in levels(cellchat@idents) ){
						try({
							#SOURCENAME="T cell" #levels(cellchat@idents)
							print(SOURCENAME)
							SOURCE = which( levels(cellchat@idents) == SOURCENAME)
							TARGET = 1:length(levels(cellchat@idents))
							#TARGET = setdiff( 1:length(unique(cellchat@idents)) ,c(SOURCE) )
							#TARGET = c(2,3,4,6,9,10)
							PLOT=netVisual_bubble(cellchat,sources.use=SOURCE,targets.use=TARGET,remove.isolate=T,font.size=12)
							OUTNAME=paste0(OUTDIR,"/","BubblePlot/",SOURCENAME,".pdf");pdf(OUTNAME,width=7,height=14,useDingbats=F);print(PLOT);dev.off()
						},silent=F)
						rm(PLOT,OUTNAME,SOURCE,TARGET,SOURCENAME)
					}
				}
				
				# circular and hierachical vizualisation 
				if(T){
					# Signif pathway
					SIGNALING=cellchat@netP$pathways
					#same: SIGNALING=unique(df.net [ df.net$pval <=0.05, "pathway_name"])
					dir.create(paste0(OUTDIR,"/pair_LR/"))
					pairLR <- extractEnrichedLR(cellchat, signaling = SIGNALING, geneLR.return = FALSE)
					for( index in 1:nrow(pairLR)){
						LR.show <- pairLR[index,] # show one ligand-receptor pair
						print(LR.show)
						#length(levels(cellchat@idents))
						vertex.receiver = seq(1,4) # a numeric vector: the cell group number TARGET
						PLOT=netVisual_individual(cellchat,signaling=SIGNALING,pairLR.use=LR.show,vertex.receiver=vertex.receiver,layout="hierarchy") #hierarchy, circle, chord
						#PLOT0netVisual_individual(cellchat,signaling=SIGNALING,pairLR.use=LR.show,vertex.receiver=vertex.receiver,layout="circle") #hierarchy, circle, chord
						#df.net[ df.net$interaction_name == LR.show, ]
						
						dev.off()
						
						OUTNAME=paste0(OUTDIR,"/pair_LR/",LR.show,".pdf");pdf(OUTNAME,width=14,height=7,useDingbats=F);print(PLOT);dev.off()
					}
					rm(PLOT,OUTNAME,vertex.receiver,LR.show,pairLR)
				}
			}else{
				print("NO SIGNIF INTERRACTIONS")
			}
			#clean
			rm(cellchat,SIGNALING)
		}
		rm(GEX,ORIGIN,SELECT,OUTDIR,df.net,PPI,DATABASE)
	}
	options(stringsAsFactors = TRUE)
}
