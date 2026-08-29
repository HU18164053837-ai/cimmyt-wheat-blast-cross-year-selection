options(stringsAsFactors=FALSE, warn=1)
set.seed(20260811)
root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
base <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
indir <- file.path(base, "history", "cross_year_analysis", "tables")
out <- file.path(base, "history", "evidence_pedigree_analysis")
tabdir <- file.path(out, "tables"); figdir <- file.path(out, "figures"); repdir <- file.path(out, "reports"); logdir <- file.path(out, "logs")
for (p in c(tabdir,figdir,repdir,logdir)) dir.create(p, recursive=TRUE, showWarnings=FALSE)

eligible <- read.csv(file.path(indir,"06_eligible_cross_year_gid_candidates.csv"), check.names=FALSE, encoding="UTF-8")
gny <- read.csv(file.path(indir,"04_gid_nursery_year_evidence.csv"), check.names=FALSE, encoding="UTF-8")
stopifnot(nrow(eligible)==104, !anyDuplicated(eligible$GID), all(eligible$n_years>=2))
write_out <- function(x,n) write.csv(x,file.path(tabdir,n),row.names=FALSE,na="",fileEncoding="UTF-8")
mode_text <- function(x) {x<-x[!is.na(x)&trimws(x)!=""]; if(!length(x)) NA_character_ else names(sort(table(x),decreasing=TRUE))[1]}

eligible$evidence_band <- ifelse(eligible$n_years==2,"2-year_provisional",ifelse(eligible$n_years<=4,"3-4-year_intermediate","5-7-year_long_term"))
eligible$decision_class <- "do_not_prioritize"
eligible$decision_class[eligible$evidence_band=="2-year_provisional" & eligible$robustness_tier=="strong"] <- "high_performance_needs_confirmation"
eligible$decision_class[eligible$evidence_band=="3-4-year_intermediate" & eligible$robustness_tier %in% c("strong","moderate")] <- "intermediate_support"
eligible$decision_class[eligible$evidence_band=="5-7-year_long_term" & eligible$robustness_tier %in% c("strong","moderate")] <- "long_term_priority"
eligible$decision_class[eligible$evidence_band=="5-7-year_long_term" & eligible$robustness_tier=="not_robust"] <- "long_term_susceptible_or_unstable_reference"
eligible <- eligible[order(factor(eligible$decision_class,levels=c("long_term_priority","intermediate_support","high_performance_needs_confirmation","long_term_susceptible_or_unstable_reference","do_not_prioritize")),eligible$median_percentile),]
write_out(eligible,"01_evidence_stratified_candidates.csv")

evidence_summary <- do.call(rbind,lapply(split(seq_len(nrow(eligible)),eligible$evidence_band),function(ii){z<-eligible[ii,]; data.frame(evidence_band=z$evidence_band[1],n_gid=nrow(z),strong=sum(z$robustness_tier=="strong"),moderate=sum(z$robustness_tier=="moderate"),not_robust=sum(z$robustness_tier=="not_robust"),median_of_median_percentile=median(z$median_percentile),median_q90_disease=median(z$q90_disease))}))
write_out(evidence_summary,"02_evidence_band_summary.csv")

longterm <- eligible[eligible$n_years>=5,]
write_out(longterm,"03_long_term_repeated_gid_assessment.csv")
priority <- eligible[eligible$decision_class %in% c("long_term_priority","intermediate_support","high_performance_needs_confirmation"),]
write_out(priority,"04_priority_validation_portfolio.csv")

# Aggregate nursery-year evidence to one GID-year trajectory point; no complete panel is constructed.
keep <- gny$GID %in% eligible$GID
gy <- do.call(rbind,lapply(split(which(keep),interaction(gny$GID[keep],gny$year[keep],drop=TRUE)),function(ii){z<-gny[ii,]; data.frame(GID=z$GID[1],year=z$year[1],n_nurseries=length(unique(z$nursery)),n_observations=sum(z$n_observations),median_disease=median(z$median_disease),q90_disease=max(z$q90_disease),median_percentile=median(z$median_percentile),q75_percentile=max(z$q75_percentile))}))
write_out(gy,"05_gid_year_trajectories.csv")

parse_ancestors <- function(cross){
  if(is.na(cross)||trimws(cross)=="") return(character())
  x <- toupper(cross)
  x <- gsub("\\([^)]*\\)","",x)
  parts <- unlist(strsplit(x,"/+"))
  parts <- trimws(parts)
  parts <- gsub("^[0-9]+\\*","",parts)
  parts <- gsub("\\*[0-9]+$","",parts)
  parts <- gsub("[^A-Z0-9.# _-]","",parts)
  parts <- gsub("[[:space:]]+"," ",parts)
  unique(parts[nchar(parts)>=2 & !grepl("^[0-9]+$",parts)])
}
sets <- setNames(lapply(eligible$Cross,parse_ancestors),as.character(eligible$GID))
token_rows <- do.call(rbind,lapply(names(sets),function(id) if(length(sets[[id]])) data.frame(GID=id,ancestor_token=sets[[id]]) else data.frame(GID=id,ancestor_token=NA_character_)))
write_out(token_rows,"06_gid_ancestor_tokens.csv")

ids <- names(sets); sim <- matrix(0,length(ids),length(ids),dimnames=list(ids,ids))
for(i in seq_along(ids)) for(j in i:length(ids)){
  a<-sets[[i]]; b<-sets[[j]]; u<-union(a,b)
  s<-if(!length(u)) 0 else length(intersect(a,b))/length(u)
  sim[i,j]<-sim[j,i]<-s
}
diag(sim)<-1
pairs <- do.call(rbind,lapply(seq_len(length(ids)-1),function(i) do.call(rbind,lapply((i+1):length(ids),function(j) data.frame(GID1=ids[i],GID2=ids[j],jaccard_similarity=sim[i,j],shared_ancestors=paste(intersect(sets[[i]],sets[[j]]),collapse=";"))))))
pairs <- pairs[order(-pairs$jaccard_similarity),]
write_out(pairs,"07_pairwise_pedigree_similarity.csv")

hc <- hclust(as.dist(1-sim),method="average")
eligible$pedigree_cluster_sim30 <- as.integer(cutree(hc,h=.70)[as.character(eligible$GID)])
eligible$pedigree_cluster_sim40 <- as.integer(cutree(hc,h=.60)[as.character(eligible$GID)])
write_out(eligible,"08_candidates_with_pedigree_clusters.csv")
sens <- data.frame(minimum_similarity=c(.20,.30,.40),n_clusters=sapply(c(.80,.70,.60),function(h)length(unique(cutree(hc,h=h)))),largest_cluster=sapply(c(.80,.70,.60),function(h)max(table(cutree(hc,h=h)))))
write_out(sens,"09_cluster_threshold_sensitivity.csv")

freq <- sort(table(unlist(sets)),decreasing=TRUE)
ancestry <- data.frame(ancestor_token=names(freq),n_gid=as.integer(freq),pct_eligible_gid=100*as.integer(freq)/nrow(eligible))
write_out(ancestry,"10_ancestor_frequency.csv")

strong_ids <- as.character(eligible$GID[eligible$robustness_tier=="strong"])
strong_pairs <- pairs[pairs$GID1 %in% strong_ids & pairs$GID2 %in% strong_ids,]
write_out(strong_pairs,"11_strong_candidate_pedigree_similarity.csv")
strong <- eligible[eligible$robustness_tier=="strong",]
strong$max_similarity_to_other_strong <- sapply(as.character(strong$GID),function(id){o<-setdiff(strong_ids,id); if(!length(o))0 else max(sim[id,o])})
strong$independent_source_flag <- strong$max_similarity_to_other_strong < .30
write_out(strong,"12_strong_candidate_independence.csv")

cols <- c("2-year_provisional"="#4F718C","3-4-year_intermediate"="#C9A227","5-7-year_long_term"="#3F6B3A")
draw_evidence <- function(){
  layout(matrix(c(1,2),1,2),widths=c(1.15,.85)); par(mar=c(3.5,3.8,2.2,.8),mgp=c(2.1,.65,0),tcl=-.25,family="sans",cex=.8)
plot(eligible$n_years,eligible$q90_disease,pch=21,bg=cols[eligible$evidence_band],col="white",cex=pmin(2.1,.75+.10*eligible$n_nurseries),xlim=c(1.75,7.35),ylim=c(-2,104),xlab="Years with observations",ylab="90th percentile disease index (%)")
  abline(h=c(10,20),lty=3,col="#777777")
title("a  Evidence duration and adverse-tail disease",adj=0,font.main=2,cex.main=.92)
legend(3.15,58,legend=c("2-year provisional","3–4-year intermediate","5–7-year long-term"),pt.bg=cols,pch=21,bty="n",cex=.62)
  bands <- c("2-year_provisional","3-4-year_intermediate","5-7-year_long_term")
  mat <- sapply(bands,function(b)table(factor(eligible$robustness_tier[eligible$evidence_band==b],levels=c("strong","moderate","not_robust"))))
barplot(mat,beside=FALSE,col=c("#3F6B3A","#C9A227","#B8B2A7"),border=NA,names.arg=c("2 y","3–4 y","5–7 y"),ylab="Number of repeated GIDs")
legend("topright",legend=rownames(mat),fill=c("#3F6B3A","#C9A227","#B8B2A7"),bty="n",cex=.60)
title("b  Robustness by evidence duration",adj=0,font.main=2,cex.main=.90)
}

focus_ids <- unique(c(strong_ids,as.character(eligible$GID[eligible$n_years>=3 & eligible$robustness_tier=="moderate"])))
focus_sim <- sim[focus_ids,focus_ids,drop=FALSE]
ord <- order.dendrogram(as.dendrogram(hclust(as.dist(1-focus_sim),method="average")))
focus_sim <- focus_sim[ord,ord,drop=FALSE]
draw_pedigree <- function(){
  layout(matrix(c(1,2),1,2),widths=c(1.25,.75)); par(family="sans",cex=.75)
par(mar=c(5.5,5.5,2.2,1)); image(seq_len(nrow(focus_sim)),seq_len(ncol(focus_sim)),t(focus_sim[nrow(focus_sim):1,]),col=colorRampPalette(c("#F7F5EC","#D7C98D","#3F6B3A"))(50),axes=FALSE,xlab="",ylab="",zlim=c(0,1))
  axis(1,at=seq_len(nrow(focus_sim)),labels=colnames(focus_sim),las=2,cex.axis=.5)
  axis(2,at=seq_len(ncol(focus_sim)),labels=rev(rownames(focus_sim)),las=2,cex.axis=.5)
  title("a  Pedigree Jaccard similarity",adj=0,font.main=2,cex.main=.95)
par(mar=c(3.0,5.5,2.2,1.4)); topa<-head(ancestry,15); barplot(rev(topa$n_gid),names.arg=rev(topa$ancestor_token),horiz=TRUE,las=1,col="#3F6B3A",border=NA,xlab="",cex.names=.55)
  title("b  Frequent pedigree components",adj=0,font.main=2,cex.main=.95)
}

# Publication-safe Arial/Helvetica/sans fallback is kept explicit for source preflight.
export_plot <- function(fun,stem,width=183,height=100){w<-width/25.4;h<-height/25.4;if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(file.path(figdir,paste0(stem,".svg")),width=w,height=h,system_fonts=list(sans="Arial"))else svg(file.path(figdir,paste0(stem,".svg")),width=w,height=h,family="sans");fun();dev.off();cairo_pdf(file.path(figdir,paste0(stem,".pdf")),width=w,height=h,family="sans");fun();dev.off();if(requireNamespace("ragg",quietly=TRUE))ragg::agg_tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw")else tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",family="sans");fun();dev.off();png(file.path(figdir,paste0(stem,".png")),width=w,height=h,units="in",res=300,family="sans");fun();dev.off()}
export_plot(draw_evidence,"Figure_3_evidence_duration_and_robustness")
export_plot(draw_pedigree,"Figure_4_pedigree_similarity_and_ancestry",height=125)

long_priority <- eligible[eligible$decision_class=="long_term_priority",]
intermediate <- eligible[eligible$decision_class=="intermediate_support",]
confirm <- eligible[eligible$decision_class=="high_performance_needs_confirmation",]
independent_n <- sum(strong$independent_source_flag)
report <- paste0("# 跨年度证据分层与候选材料系谱来源分析\n\n",
"## 核心判断\n\n",
"104个真实重复GID中，",sum(eligible$n_years==2),"个仅有2年证据，",sum(eligible$n_years%in%3:4),"个有3–4年证据，",sum(eligible$n_years>=5),"个有5–7年证据。上一阶段8个strong候选全部属于2年初步证据，因此定位为高表现、待确认材料，而不是长期稳定材料。\n\n",
"长期证据层中，",nrow(long_priority),"个材料同时保持moderate或更高稳健性：",paste(long_priority$GID,collapse="、"),"。中期证据层有",nrow(intermediate),"个材料可进入扩大验证：",paste(intermediate$GID,collapse="、"),"。\n\n",
"## 系谱来源\n\n",
"系谱被拆分为标准化亲本/材料名称集合，并用Jaccard相似度描述共同组成。相似度≥0.30仅作为探索性来源簇阈值；阈值敏感性见`09_cluster_threshold_sensitivity.csv`。8个strong候选中，",independent_n,"个与其他strong候选的最大系谱相似度低于0.30，可视为相对独立的候选来源；其余材料可能共享较多育种背景。\n\n",
"## 决策建议\n\n",
"1. 将长期优先材料作为跨年稳定性骨架；\n2. 将8个2年strong材料安排到相同高区分力环境中复验；\n3. 在复验组合中优先保留系谱相似度低的材料，避免验证资源集中在同一背景；\n4. 对长期感病或不稳定材料保留为环境反应参照，不作为抗性候选。\n\n",
"## 限制\n\n系谱字符串不等同于标记数据或真实亲缘系数；名称拆分、别名和复杂回交表达会影响Jaccard相似度。因此聚类用于候选组合设计，不能替代基因型亲缘分析。\n")
writeLines(report,file.path(repdir,"evidence_pedigree_report.md"),useBytes=TRUE)
qa <- "# QA\n\n- 输入104个跨年度合格GID，无额外删除。\n- 轨迹表按GID×年份汇总，不构造缺失年份或完整面板。\n- 系谱聚类使用全部104个GID；图中热图聚焦8个strong及有≥3年证据的moderate材料。\n- 所有图件由R后端生成，宽度183 mm，含SVG/PDF/600 dpi TIFF/PNG。\n"
writeLines(qa,file.path(repdir,"QA_notes.md"),useBytes=TRUE)
sink(file.path(logdir,"sessionInfo.txt"));print(sessionInfo());sink()
summary <- data.frame(metric=c("eligible_gid","two_year_gid","three_four_year_gid","five_seven_year_gid","high_performance_needs_confirmation","intermediate_support","long_term_priority","independent_strong_sources"),value=c(nrow(eligible),sum(eligible$n_years==2),sum(eligible$n_years%in%3:4),sum(eligible$n_years>=5),nrow(confirm),nrow(intermediate),nrow(long_priority),independent_n))
write_out(summary,"00_analysis_summary.csv");cat(paste(summary$metric,summary$value,sep="="),sep="\n")
