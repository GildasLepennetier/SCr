FAM<-function(){
	for ( NAMEIDENT in NAMES_CLU ){ 
		print(paste("DOING:",NAMEIDENT))
		Idents(GEX)<-NAMEIDENT
		allMarkers <- FindAllMarkers(GEX)
		write.xlsx(allMarkers, paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".xlsx"))
		saveRDS(allMarkers,paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
		avgexp <- AverageExpression( GEX[allMarkers$gene] ) ; avgexp$RNA$gene=rownames(avgexp$RNA) ; avgexp=avgexp$RNA
		allMarkersExpr=merge.data.frame(x=avgexp,y=allMarkers,by="gene") %>% arrange(cluster,p_val_adj )
		write.xlsx(as.data.frame(allMarkersExpr %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=20)),paste0(SUBDIR,"/allMarkers+avgexp+top20.",NAMEIDENT,".xlsx"))
	}
	rm(allMarkers,avgexp,allMarkersExpr)
}
