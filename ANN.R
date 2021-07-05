ANN<-function(){
	for ( NAMEIDENT in NAMES_CLU ) { # c("RNA_snn_res.0.1", "RNA_snn_res.0.2") ) { #
		Idents(GEX) <- NAMEIDENT
		print(paste0("doing: ",NAMEIDENT, " (",length(unique(Idents(GEX)))," clusters)"))
		# # # # CIPR - Cluster Identity PRedictor
		if(T){
			#required for CIPR
			allMarkers<-readRDS(paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
			
			DATABASE="mmrnaseq"
			#DATABASE="immgen"
			
			# Plot summarizing top scoring references per cluster (logFC comparison)
			PLOT <- CIPR(input_dat = allMarkers,comp_method="logfc_dot_product",reference=DATABASE,top_num=1,keep_top_var=50) + 
				grids(axis="x",color="grey92",size=NULL,linetype=NULL) #"dotted")
			#OUTNAME=paste0("7.CIPR.",NAMEIDENT,"_",DATABASE,".pdf");pdf(OUTNAME,width=14,height=7,useDingbats=F);print(PLOT);dev.off()
			dev.off()

			ANNOT = CIPR_top_results %>% 
				select( cluster,reference_cell_type) %>% 
				distinct() %>%
				group_by(cluster) %>% mutate( annot = paste0(reference_cell_type,collapse="+")) %>%
				select(cluster,annot) %>% distinct()
			#Deal with unannoted clusters
			NOT_ANNOT = levels(Idents(GEX)) [ ! levels(Idents(GEX)) %in% ANNOT$cluster ]
			if(length(NOT_ANNOT)>0){
				tmp_df=data.frame(cluster=NOT_ANNOT, annot="unknown")
				ANNOT=rbind(ANNOT,tmp_df) ; rm(tmp_df)
			}
			#USE THE ROWNAME TO AVOID ERROR IN CLUSTER NAME ~ LEVELS
			ANNOT=as.data.frame(ANNOT[order(ANNOT$cluster),])
			rownames(ANNOT)=as.character(ANNOT$cluster)
			write.xlsx(ANNOT, paste0(SUBDIR,"/CIPR.",NAMEIDENT,".xlsx"))
			#ADD in the data
			NEWNAME=paste0("annot_",NAMEIDENT,"_CIPR")
			GEX@meta.data[,NEWNAME] <- GEX@meta.data[,NAMEIDENT]
			levels(GEX@meta.data[,NEWNAME]) <- ANNOT[ levels(GEX@meta.data[,NEWNAME]),"annot"]
			rm(CIPR_top_results,CIPR_all_results, top_plots,PLOT,allMarkers,ANNOT)
			rm(DATABASE,NEWNAME,NOT_ANNOT)
			
		} 
	}
}
	
