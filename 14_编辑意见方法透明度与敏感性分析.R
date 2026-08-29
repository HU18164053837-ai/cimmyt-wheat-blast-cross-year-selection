options(stringsAsFactors=FALSE)
root <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
out <- file.path(root,"history","editor_requested_transparency_analysis")
dir.create(file.path(out,"tables"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"figures"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"reports"),recursive=TRUE,showWarnings=FALSE)

env <- read.csv(file.path(root,"history","cross_year_analysis","tables","02_environment_pressure_and_discrimination.csv"),check.names=FALSE)
rob <- read.csv(file.path(root,"history","cross_year_analysis","tables","05_gid_cross_year_robustness.csv"),check.names=FALSE)
elig <- read.csv(file.path(root,"history","cross_year_analysis","tables","06_eligible_cross_year_gid_candidates.csv"),check.names=FALSE)
priority <- read.csv(file.path(root,"history","evidence_pedigree_analysis","tables","04_priority_validation_portfolio.csv"),check.names=FALSE)
portfolio <- read.csv(file.path(root,"history","validation_portfolio_analysis","tables","01_nonredundant_validation_portfolio.csv"),check.names=FALSE)

prank <- function(x) (rank(x,ties.method="average")-1)/max(1,length(x)-1)
env$pressure_component <- prank(env$balanced_mean)
env$discrimination_component <- prank(env$discrimination_iqr)
env$coverage_component <- prank(env$n_nurseries)
env$saturation_pct <- env$floor_pct_balanced + env$ceiling_pct_balanced
env$nonsaturation_component <- 1-prank(env$saturation_pct)
cn <- c("pressure_component","discrimination_component","coverage_component","nonsaturation_component")
env$equal_weight_score <- rowMeans(env[,cn])
env$equal_weight_rank <- rank(-env$equal_weight_score,ties.method="min")

# Leave-one-component-out and prespecified alternative weights.
variants <- list(
  equal=c(.25,.25,.25,.25),
  no_pressure=c(0,1/3,1/3,1/3),
  no_discrimination=c(1/3,0,1/3,1/3),
  no_coverage=c(1/3,1/3,0,1/3),
  no_nonsaturation=c(1/3,1/3,1/3,0),
  discrimination_emphasis=c(.20,.40,.20,.20),
  pressure_emphasis=c(.40,.20,.20,.20)
)
sens <- do.call(rbind,lapply(names(variants),function(nm){
  w <- variants[[nm]]; sc <- as.matrix(env[,cn])%*%w; rk <- rank(-sc,ties.method="min")
  data.frame(variant=nm,pressure_weight=w[1],discrimination_weight=w[2],coverage_weight=w[3],nonsaturation_weight=w[4],
             spearman_with_equal=cor(sc,env$equal_weight_score,method="spearman"),
             top6_overlap=sum(order(-sc)[1:6] %in% order(-env$equal_weight_score)[1:6]),
             top10_overlap=sum(order(-sc)[1:10] %in% order(-env$equal_weight_score)[1:10]))
}))
write.csv(sens,file.path(out,"tables","01_environment_score_sensitivity.csv"),row.names=FALSE)

rankings <- do.call(rbind,lapply(names(variants),function(nm){
  w<-variants[[nm]]; sc<-drop(as.matrix(env[,cn])%*%w)
  data.frame(variant=nm,environment=env$environment,score=sc,rank=rank(-sc,ties.method="min"))
}))
write.csv(rankings,file.path(out,"tables","02_environment_rankings_all_variants.csv"),row.names=FALSE)

# Equal-nursery versus pooled summaries.
env$balanced_rank <- rank(-env$balanced_mean,ties.method="average")
env$pooled_rank <- rank(-env$pooled_mean,ties.method="average")
env$rank_shift <- env$balanced_rank-env$pooled_rank
cmp <- data.frame(
  n_environments=nrow(env),
  spearman_balanced_vs_pooled=cor(env$balanced_mean,env$pooled_mean,method="spearman"),
  median_absolute_mean_difference=median(abs(env$balanced_mean-env$pooled_mean)),
  maximum_absolute_mean_difference=max(abs(env$balanced_mean-env$pooled_mean)),
  n_rank_shift_ge_3=sum(abs(env$rank_shift)>=3),
  n_rank_shift_ge_5=sum(abs(env$rank_shift)>=5)
)
write.csv(cmp,file.path(out,"tables","03_balanced_vs_pooled_summary.csv"),row.names=FALSE)
write.csv(env[order(-abs(env$rank_shift)),c("environment","n_nurseries","balanced_mean","pooled_mean","balanced_minus_pooled_mean","balanced_rank","pooled_rank","rank_shift")],file.path(out,"tables","04_balanced_vs_pooled_environments.csv"),row.names=FALSE)

# Auditable selection funnel from existing locked tables.
funnel <- data.frame(
  stage=c("Linkable GIDs in primary long table","Cross-year eligibility (>=2 years, >=2 nurseries, >=4 observations)","Strong or moderate robustness tier","Priority evidence set before pedigree de-duplication","Nonredundant prospective resistance portfolio"),
  n=c(4088,nrow(elig),sum(elig$robustness_tier %in% c("strong","moderate")),nrow(priority),sum(portfolio$selected_for_resistance_validation)),
  source=c("Primary audit table","06_eligible_cross_year_gid_candidates.csv","06_eligible_cross_year_gid_candidates.csv","04_priority_validation_portfolio.csv","01_nonredundant_validation_portfolio.csv")
)
write.csv(funnel,file.path(out,"tables","05_candidate_selection_funnel.csv"),row.names=FALSE)

# Supplementary figure: reviewer-facing transparency checks.
cols <- c(green="#3F6B3A",gold="#C9A227",brown="#8B5E3C",blue="#4F718C",grey="#B8B2A7")
ss <- sens[sens$variant!="equal",]
draw <- function(){
  layout(matrix(1:3,1,3),widths=c(1,1,1.05)); par(family="sans",cex=.78)
  par(mar=c(4.2,4.2,2.6,.7)); plot(env$pooled_mean,env$balanced_mean,pch=21,bg=cols["green"],col="white",xlab="Pooled disease mean (%)",ylab="Equal-nursery mean (%)"); abline(0,1,lty=2,col=cols["brown"],lwd=1.5); title("a  Weighting comparison",adj=0,font.main=2,cex.main=.95)
  par(mar=c(7,4.2,2.6,.7)); bp<-barplot(ss$spearman_with_equal,names.arg=gsub("_","\n",ss$variant),las=2,ylim=c(0,1),col=c(cols["brown"],cols["gold"],cols["blue"],cols["grey"],cols["green"],cols["gold"]),border=NA,ylab="Spearman correlation"); text(bp,ss$spearman_with_equal,sprintf("%.2f",ss$spearman_with_equal),pos=3,cex=.62); title("b  Score sensitivity",adj=0,font.main=2,cex.main=.95)
  par(mar=c(6.4,4.2,2.6,.7)); barplot(funnel$n,names.arg=c("Linkable","Eligible","Robust","Priority","Final"),las=2,log="y",col=c(cols["blue"],cols["grey"],cols["gold"],cols["brown"],cols["green"]),border=NA,ylab="Number of GIDs (log scale)"); text(seq(.7,by=1.2,length.out=nrow(funnel)),funnel$n,labels=funnel$n,pos=3,cex=.65); title("c  Candidate decision funnel",adj=0,font.main=2,cex.main=.95)
}
base<-file.path(out,"figures","Figure_S4_method_transparency_checks")
png(paste0(base,".png"),width=3600,height=1450,res=400);draw();dev.off()
tiff(paste0(base,".tiff"),width=7.2,height=2.9,units="in",res=600,compression="lzw");draw();dev.off()
if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(paste0(base,".svg"),width=7.2,height=2.9,system_fonts=list(sans="Arial")) else svg(paste0(base,".svg"),width=7.2,height=2.9,family="sans");draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=7.2,height=2.9,family="sans");draw();dev.off()

writeLines(c("# Editor-requested method transparency checks",sprintf("Equal-nursery and pooled disease means had Spearman rho %.3f; %d of %d environments shifted by at least three rank positions.",cmp$spearman_balanced_vs_pooled,cmp$n_rank_shift_ge_3,cmp$n_environments),sprintf("Across six alternative component specifications, correlation with the equal-weight information score ranged from %.3f to %.3f and top-six overlap ranged from %d to %d environments.",min(ss$spearman_with_equal),max(ss$spearman_with_equal),min(ss$top6_overlap),max(ss$top6_overlap)),paste("Candidate funnel:",paste(paste(funnel$n,funnel$stage),collapse="; "))),file.path(out,"reports","analysis_report.md"))
capture.output(sessionInfo(),file=file.path(out,"reports","sessionInfo.txt"))
print(cmp);print(sens);print(funnel)
