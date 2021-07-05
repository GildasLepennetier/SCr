MEXL<-function(){
	for ( NAMEIDENT in NAMES_CLU ){
		allMarkersExpr <- readRDS(paste0(SUBDIR,"/allMarkers.",NAMEIDENT,".Rds"))
		Annotation <-read.xlsx(paste0(SUBDIR,"/CIPR.",NAMEIDENT,".xlsx"))
		MERGED=merge.data.frame(allMarkersExpr,Annotation,by="cluster")
		write.xlsx(MERGED,paste0(SUBDIR,"/allMarkers.",NAMEIDENT,"_CIPR.xlsx"))
		rm(allMarkersExpr,Annotation,MERGED)
	}
}
