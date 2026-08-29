base <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
out <- file.path(base,"history","temporal_uncertainty_negative_control","figures","Figure_research_workflow")
cols <- c("#EEF1E7", "#DCE7D5", "#F3E6BE", "#E8D8C8", "#DCE5E8")
publication_font_family <- "Times New Roman"
draw <- function() {
  par(mar=c(.2,.2,.2,.2), family="serif"); plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
  boxes <- list(c(.02,.69,.19,.91),c(.22,.69,.39,.91),c(.42,.69,.59,.91),c(.62,.69,.79,.91),c(.82,.69,.98,.91),
                c(.12,.25,.34,.50),c(.39,.25,.61,.50),c(.66,.25,.88,.50))
  labels <- c("Public nursery archive\n2018–2024","Compatibility audit\nand duplicate control",
              "Within-environment\npercentile ranks","Cross-year evidence\nwithout panel completion",
              "14-candidate\nnonredundant portfolio","Common-check and\nbackground adjustment",
              "Past-only next-year\nenrichment tests","Prospective balanced\nvalidation design")
  sub <- c("24 phenotype files","15HZAN excluded","structural absence\nretained","GID-linked recurrence",
           "breeding roles assigned","environment-cluster\nuncertainty","cluster bootstrap\n+ permutation","recommended next step")
  for(i in seq_along(boxes)) {
    b <- boxes[[i]]; polygon(c(b[1],b[3],b[3],b[1]),c(b[2],b[2],b[4],b[4]),col=cols[(i-1)%%5+1],border="#3F5145",lwd=1.2)
    text(mean(b[c(1,3)]),mean(b[c(2,4)])+.035,labels[i],font=2,cex=.62)
    text(mean(b[c(1,3)]),b[2]+.038,sub[i],cex=.50,col="#4D5559")
  }
  arrow <- function(x0,y0,x1,y1) arrows(x0,y0,x1,y1,length=.08,lwd=1.4,col="#3F5145")
  for(i in 1:4) arrow(boxes[[i]][3],.80,boxes[[i+1]][1],.80)
  segments(.90,.69,.90,.58,lwd=1.4,col="#3F5145"); segments(.90,.58,.23,.58,lwd=1.4,col="#3F5145"); arrow(.23,.58,.23,.50)
  arrow(.34,.375,.39,.375); arrow(.61,.375,.66,.375)
  text(.53,.625,"Retrospective evidence integration",cex=.72,col="#3F5145")
  text(.5,.10,"Scope: candidate enrichment and risk triage—not independent external validation",font=3,cex=.72,col="#8B5E3C")
}
png(paste0(out,".png"),width=3600,height=1900,res=400); draw(); dev.off()
tiff(paste0(out,".tiff"),width=7.2,height=3.8,units="in",res=600,compression="lzw"); draw(); dev.off()
if (requireNamespace("svglite", quietly=TRUE)) svglite::svglite(paste0(out,".svg"),width=7.2,height=3.8,system_fonts=list(serif="Times New Roman")) else svg(paste0(out,".svg"),width=7.2,height=3.8,family="serif",pointsize=8); draw(); dev.off()
cairo_pdf(paste0(out,".pdf"),width=7.2,height=3.8,family="serif",pointsize=8); draw(); dev.off()
