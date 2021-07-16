DOT<-function(){
	FEATURES<-c(      "Cd8a","Cd4","Col1a1","Krt18","Pecam1", "S100a8",  "Cd68", "Cx3cr1",  "Olig1",  "Gja1")
	names(FEATURES)=c("CD8", "CD4","Fibro.","Epit.","Endo.",   "Granulo.","Mono.","Micro.","Oligo.","Astro.")

	NAMEIDENT = "CellTypeFinal"
	#DefaultAssay(GEX)
	Idents(GEX)<-NAMEIDENT
	IDENTS = unique(Idents(GEX))
	PLOT=DotPlot(GEX,features=FEATURES,assay="RNA",idents=IDENTS,dot.scale=10) + 
		theme_minimal(base_size=20) + 
		theme(axis.text.x=element_text(angle=90)) + 
		ggtitle("Gene markers")
	OUTNAME=paste0(SUBDIR,"/DotPlot_GeneMarkers.pdf") ; pdf(OUTNAME,width=20,height=8) ; print(PLOT);dev.off()

}
