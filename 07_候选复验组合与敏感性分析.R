options(stringsAsFactors=FALSE,warn=1)
set.seed(20260811)
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);base<-Sys.getenv("CIMMYT_ANALYSIS_ROOT",unset=root)
longfile<-file.path(base,"history","standardized","historical_long_primary_observed.csv")
crossdir<-file.path(base,"history","cross_year_analysis","tables")
pedir<-file.path(base,"history","evidence_pedigree_analysis","tables")
out<-file.path(base,"history","validation_portfolio_analysis");tabdir<-file.path(out,"tables");figdir<-file.path(out,"figures");repdir<-file.path(out,"reports");logdir<-file.path(out,"logs")
for(p in c(tabdir,figdir,repdir,logdir))dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_out<-function(x,n)write.csv(x,file.path(tabdir,n),row.names=FALSE,na="",fileEncoding="UTF-8")
qfun<-function(x,p)unname(quantile(x,p,type=7,na.rm=TRUE));prank<-function(x){if(length(unique(x))==1)rep(.5,length(x))else(rank(x,ties.method="average")-1)/(length(x)-1)}

d<-read.csv(longfile,check.names=FALSE,na.strings=c("","NA"),encoding="UTF-8-BOM")
env<-read.csv(file.path(crossdir,"02_environment_pressure_and_discrimination.csv"),check.names=FALSE)
cand<-read.csv(file.path(pedir,"08_candidates_with_pedigree_clusters.csv"),check.names=FALSE)
stopifnot(!any(d$nursery=="15HZAN"),nrow(d)==30475,nrow(env)==38,nrow(cand)==104)
d$GID<-trimws(as.character(d$GID));d$GID[d$GID==""]<-NA_character_;d$year<-as.integer(d$year);d$disease_index_pct<-as.numeric(d$disease_index_pct)

# Build a non-redundant testing portfolio. Within each strong pedigree cluster retain the better adverse-tail profile.
anchors<-cand[cand$decision_class=="long_term_priority",]
intermediate<-cand[cand$decision_class=="intermediate_support",]
strong<-cand[cand$decision_class=="high_performance_needs_confirmation",]
strong<-strong[order(strong$pedigree_cluster_sim30,strong$q90_disease,strong$worst_disease,strong$median_percentile),]
strong_rep<-strong[!duplicated(strong$pedigree_cluster_sim30),]
controls<-cand[as.character(cand$GID)%in%c("304660","6387","911521"),]
make_role<-function(z,role){z$portfolio_role<-role;z}
portfolio<-rbind(make_role(anchors,"long_term_anchor"),make_role(intermediate,"intermediate_candidate"),make_role(strong_rep,"provisional_strong_unique_pedigree"),make_role(controls,"susceptible_or_variable_reference"))
portfolio$selected_for_resistance_validation<-portfolio$portfolio_role!="susceptible_or_variable_reference"
portfolio<-portfolio[order(factor(portfolio$portfolio_role,levels=c("long_term_anchor","intermediate_candidate","provisional_strong_unique_pedigree","susceptible_or_variable_reference")),portfolio$median_percentile),]
write_out(portfolio,"01_nonredundant_validation_portfolio.csv")

# Rank historical environments for informativeness, not for expected future performance.
env$pressure_component<-prank(env$balanced_mean)
env$discrimination_component<-prank(env$discrimination_iqr)
env$coverage_component<-prank(env$n_nurseries)
env$saturation_pct<-env$floor_pct_balanced+env$ceiling_pct_balanced
env$nonsaturation_component<-1-prank(env$saturation_pct)
env$historical_information_score<-rowMeans(env[,c("pressure_component","discrimination_component","coverage_component","nonsaturation_component")])
env$information_rank<-rank(-env$historical_information_score,ties.method="min")
env<-env[order(env$information_rank),]
write_out(env,"02_historical_environment_information_ranking.csv")

stratum<-do.call(rbind,lapply(split(seq_len(nrow(env)),paste(env$location,env$sowing,sep="|")),function(ii){z<-env[ii,];data.frame(location=z$location[1],sowing=z$sowing[1],n_years=nrow(z),median_balanced_mean=median(z$balanced_mean),median_balanced_median=median(z$balanced_median),median_discrimination_iqr=median(z$discrimination_iqr),median_saturation_pct=median(z$saturation_pct),median_information_score=median(z$historical_information_score),best_historical_environment=z$environment[which.max(z$historical_information_score)])}))
stratum<-stratum[order(-stratum$median_information_score),];write_out(stratum,"03_location_sowing_validation_value.csv")

# Candidate responses remain observed-only. Percentiles are calculated within nursery x environment.
g_ne<-interaction(d$nursery,d$environment,drop=TRUE)
d$within_environment_percentile<-ave(d$disease_index_pct,g_ne,FUN=function(x)rank(x,ties.method="average")/length(x))
pd<-d[!is.na(d$GID)&d$GID%in%as.character(portfolio$GID),]
g_ps<-interaction(pd$GID,pd$location,pd$sowing,drop=TRUE)
profile<-do.call(rbind,lapply(split(seq_len(nrow(pd)),g_ps),function(ii){z<-pd[ii,];data.frame(GID=z$GID[1],location=z$location[1],sowing=z$sowing[1],n_observations=nrow(z),n_years=length(unique(z$year)),n_nurseries=length(unique(z$nursery)),median_disease=median(z$disease_index_pct),q90_disease=qfun(z$disease_index_pct,.90),worst_disease=max(z$disease_index_pct),median_percentile=median(z$within_environment_percentile),q75_percentile=qfun(z$within_environment_percentile,.75))}))
profile<-merge(profile,portfolio[,c("GID","portfolio_role")],by="GID",all.x=TRUE)
profile$vulnerability_flag<-profile$q75_percentile>.60|profile$q90_disease>20
write_out(profile,"04_portfolio_location_sowing_profiles.csv")
write_out(profile[profile$vulnerability_flag,],"05_candidate_vulnerability_flags.csv")

classify<-function(z){
  years<-length(unique(z$year));nurs<-length(unique(z$nursery));nobs<-nrow(z)
  if(years<2||nurs<2||nobs<4)return("insufficient")
  ym<-tapply(z$within_environment_percentile,z$year,median)
  medp<-median(z$within_environment_percentile);q75p<-qfun(z$within_environment_percentile,.75);q90<-qfun(z$disease_index_pct,.90);worst<-max(z$disease_index_pct);maxym<-max(ym)
  if(medp<=.40&&q75p<=.45&&maxym<=.50&&q90<=10&&worst<=20)return("strong")
  if(medp<=.45&&q75p<=.60&&maxym<=.65&&q90<=20&&worst<=50)return("moderate")
  "not_robust"
}

# Leave-one-year-out stability is defined only for GIDs with at least three observed years.
lo_ids<-as.character(portfolio$GID[portfolio$n_years>=3&portfolio$selected_for_resistance_validation])
loyo<-do.call(rbind,lapply(lo_ids,function(id){z<-d[!is.na(d$GID)&d$GID==id,];do.call(rbind,lapply(sort(unique(z$year)),function(y){zz<-z[z$year!=y,];data.frame(GID=id,excluded_year=y,remaining_years=length(unique(zz$year)),remaining_nurseries=length(unique(zz$nursery)),remaining_observations=nrow(zz),leave_one_year_tier=classify(zz))}))}))
loyo<-merge(loyo,portfolio[,c("GID","portfolio_role","robustness_tier")],by="GID",all.x=TRUE)
write_out(loyo,"06_leave_one_year_out_results.csv")
lo_summary<-do.call(rbind,lapply(split(seq_len(nrow(loyo)),loyo$GID),function(ii){z<-loyo[ii,];data.frame(GID=z$GID[1],portfolio_role=z$portfolio_role[1],original_tier=z$robustness_tier[1],n_exclusions=nrow(z),same_tier_fraction=mean(z$leave_one_year_tier==z$robustness_tier[1]),at_least_moderate_fraction=mean(z$leave_one_year_tier%in%c("strong","moderate")),insufficient_fraction=mean(z$leave_one_year_tier=="insufficient"),tiers_observed=paste(sort(unique(z$leave_one_year_tier)),collapse=";"))}))
write_out(lo_summary,"07_leave_one_year_out_summary.csv")

# Pre-specified conservative/reference/liberal threshold sensitivity on all 104 repeated GIDs.
scenarios<-data.frame(scenario=c("conservative","reference","liberal"),median_pct=c(.35,.40,.45),q75_pct=c(.40,.45,.50),max_year_pct=c(.45,.50,.60),q90_disease=c(5,10,15),worst_disease=c(10,20,30))
scenario_counts<-do.call(rbind,lapply(seq_len(nrow(scenarios)),function(i){s<-scenarios[i,];pass<-cand$median_percentile<=s$median_pct&cand$q75_percentile<=s$q75_pct&cand$max_year_median_percentile<=s$max_year_pct&cand$q90_disease<=s$q90_disease&cand$worst_disease<=s$worst_disease;data.frame(s,selected_gid=sum(pass),selected_two_year=sum(pass&cand$n_years==2),selected_three_plus=sum(pass&cand$n_years>=3))}))
write_out(scenario_counts,"08_threshold_sensitivity.csv")

loccols<-c(Jashore="#3F6B3A",Quirusillas="#C9A227",Okinawa="#8B5E3C")
draw_environment<-function(){par(mar=c(3.7,4,2.2,1),mgp=c(2.2,.7,0),tcl=-.25,family="sans",cex=.8);plot(env$balanced_mean,env$discrimination_iqr,pch=ifelse(env$sowing=="1st",21,24),bg=loccols[env$location],col="white",cex=.8+2*env$historical_information_score,xlab="Nursery-balanced mean disease (%)",ylab="Median within-nursery IQR (%)",xlim=c(0,25),ylim=c(0,33));top<-env$information_rank<=6;short<-paste0(substr(env$location[top],1,1),env$year[top],"-",substr(env$sowing[top],1,1));text(env$balanced_mean[top],env$discrimination_iqr[top],labels=short,pos=ifelse(env$discrimination_iqr[top]>25,2,4),cex=.58);legend("topleft",legend=names(loccols),pt.bg=loccols,pch=21,bty="n",cex=.65);legend("bottomright",legend=c("1st sowing","2nd sowing"),pch=c(21,24),pt.bg="#888888",col="white",bty="n",cex=.65);title("Historically informative validation environments",adj=0,font.main=2,cex.main=.95)}

res_port<-portfolio[portfolio$selected_for_resistance_validation,]
row_ids<-as.character(res_port$GID);strata<-as.vector(outer(c("Jashore","Quirusillas","Okinawa"),c("1st","2nd"),paste,sep="|"))
hm<-matrix(NA_real_,length(row_ids),length(strata),dimnames=list(row_ids,strata))
for(i in seq_len(nrow(profile))){r<-match(as.character(profile$GID[i]),row_ids);c<-match(paste(profile$location[i],profile$sowing[i],sep="|"),strata);if(!is.na(r)&&!is.na(c))hm[r,c]<-profile$median_percentile[i]}
draw_candidate<-function(){par(mar=c(7,4.8,2.2,1),family="sans",cex=.78);image(seq_len(ncol(hm)),seq_len(nrow(hm)),t(hm[nrow(hm):1,,drop=FALSE]),col=colorRampPalette(c("#2166AC","white","#B2182B"))(80),axes=FALSE,zlim=c(0,1),xlab="",ylab="");axis(1,at=seq_len(ncol(hm)),labels=gsub("\\|","\n",colnames(hm)),las=2,cex.axis=.65);axis(2,at=seq_len(nrow(hm)),labels=rev(rownames(hm)),las=2,cex.axis=.65);title("Candidate vulnerability across location × sowing strata",adj=0,font.main=2,cex.main=.95);mtext("Blue = lower disease rank; red = higher disease rank",side=1,line=5.2,cex=.65)}

draw_loyo<-function(){z<-lo_summary[order(lo_summary$at_least_moderate_fraction),];par(mar=c(4.3,5.2,2.2,1),family="sans",cex=.8);barplot(z$at_least_moderate_fraction,names.arg=z$GID,horiz=TRUE,las=1,xlim=c(0,1),col=ifelse(z$at_least_moderate_fraction==1,"#228833","#CCBB44"),border=NA,xlab="Fraction ≥ moderate after one-year exclusion",cex.names=.65);abline(v=.8,lty=3,col="#666666");title("Leave-one-year-out robustness",adj=0,font.main=2,cex.main=.95)}

# Arial/Helvetica/sans fallback is declared; all devices use the R backend.
export_plot<-function(fun,stem,width=183,height=105){w<-width/25.4;h<-height/25.4;if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(file.path(figdir,paste0(stem,".svg")),width=w,height=h,system_fonts=list(sans="Arial"))else svg(file.path(figdir,paste0(stem,".svg")),width=w,height=h,family="sans");fun();dev.off();cairo_pdf(file.path(figdir,paste0(stem,".pdf")),width=w,height=h,family="sans");fun();dev.off();if(requireNamespace("ragg",quietly=TRUE))ragg::agg_tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw")else tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",family="sans");fun();dev.off();png(file.path(figdir,paste0(stem,".png")),width=w,height=h,units="in",res=300,family="sans");fun();dev.off()}
export_plot(draw_environment,"Figure_5_validation_environment_information",height=115)
export_plot(draw_candidate,"Figure_6_candidate_location_sowing_vulnerability",height=145)
export_plot(draw_loyo,"Figure_7_leave_one_year_out_robustness",height=115)

top_env<-head(env,6);unstable<-lo_summary[lo_summary$at_least_moderate_fraction<.8,]
report<-paste0("# 候选复验组合、环境匹配与敏感性分析\n\n",
"## 非冗余复验组合\n\n复验组合包含",sum(portfolio$selected_for_resistance_validation),"个抗性候选和",sum(!portfolio$selected_for_resistance_validation),"个感病/变异参照。抗性候选由2个长期骨架、5个中期支持材料和",nrow(strong_rep),"个非冗余短期strong系谱代表构成。系谱完全相同的7627645/7627655只保留表型尾部更优的7627645作为主代表。\n\n",
"## 历史环境信息量\n\n信息评分同时考虑病害平均压力、育种圃内IQR、育种圃覆盖数和非饱和程度。排名前6的历史环境为：",paste(top_env$environment,collapse="、"),"。该评分用于选择复验条件，不表示这些具体年份可以复制。\n\n",
"## 候选薄弱环境与逐年敏感性\n\n`04_portfolio_location_sowing_profiles.csv`按真实观测总结候选在地点×播期中的疾病尾部和相对百分位；未观测组合保持缺失。逐年剔除分析覆盖",nrow(lo_summary),"个具有至少3年证据的抗性候选，其中",nrow(unstable),"个在剔除某一年后保持moderate或更高的比例低于80%：",paste(unstable$GID,collapse="、"),"。\n\n",
"## 解释边界\n\n环境信息评分是描述性多指标排序；候选热图不是平衡试验，逐年剔除也不产生新的独立重复。建议把结果用于下一轮试验分层和资源配置，而不是替代前瞻性验证。\n")
writeLines(report,file.path(repdir,"validation_portfolio_report.md"),useBytes=TRUE)
qa<-paste0("# QA\n\n- 输入主长表30475条，确认不存在15HZAN。\n- 复验候选仅来自前一阶段证据组合；同系谱strong候选按预设尾部病害排序去冗余。\n- 候选×地点×播期仅汇总真实观测，未填补缺失单元。\n- 环境评分四个分量等权，单项均为样本内百分位；敏感性场景预先固定为conservative/reference/liberal。\n- 图件由R生成，183 mm宽，SVG/PDF/600 dpi TIFF/PNG。\n")
writeLines(qa,file.path(repdir,"QA_notes.md"),useBytes=TRUE);sink(file.path(logdir,"sessionInfo.txt"));print(sessionInfo());sink()
summary<-data.frame(metric=c("portfolio_total","resistance_candidates","reference_controls","long_term_anchors","intermediate_candidates","unique_provisional_strong","loyo_candidates","loyo_below_80pct"),value=c(nrow(portfolio),sum(portfolio$selected_for_resistance_validation),sum(!portfolio$selected_for_resistance_validation),nrow(anchors),nrow(intermediate),nrow(strong_rep),nrow(lo_summary),nrow(unstable)))
write_out(summary,"00_analysis_summary.csv");cat(paste(summary$metric,summary$value,sep="="),sep="\n")
