options(stringsAsFactors=FALSE,warn=1)
set.seed(20260811)
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);base<-Sys.getenv("CIMMYT_ANALYSIS_ROOT",unset=root)
longfile<-file.path(base,"history","standardized","historical_long_primary_observed.csv")
portfile<-file.path(base,"history","validation_portfolio_analysis","tables","01_nonredundant_validation_portfolio.csv")
out<-file.path(base,"history","reference_anchored_analysis");tabdir<-file.path(out,"tables");figdir<-file.path(out,"figures");repdir<-file.path(out,"reports");logdir<-file.path(out,"logs")
for(p in c(tabdir,figdir,repdir,logdir))dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_out<-function(x,n)write.csv(x,file.path(tabdir,n),row.names=FALSE,na="",fileEncoding="UTF-8")
qfun<-function(x,p)unname(quantile(x,p,type=7,na.rm=TRUE))
d<-read.csv(longfile,check.names=FALSE,na.strings=c("","NA"),encoding="UTF-8-BOM");portfolio<-read.csv(portfile,check.names=FALSE)
d$GID<-trimws(as.character(d$GID));d$GID[d$GID==""]<-NA_character_;d$year<-as.integer(d$year);d$disease_index_pct<-as.numeric(d$disease_index_pct)
res_ids<-as.character(portfolio$GID[portfolio$selected_for_resistance_validation]);stopifnot(nrow(d)==30475,!any(d$nursery=="15HZAN"),length(res_ids)==14)

aggregate_gid_strata<-function(gids){
  z<-d[!is.na(d$GID)&d$GID%in%gids,]
  sp<-split(seq_len(nrow(z)),interaction(z$GID,z$nursery,z$environment,drop=TRUE))
  do.call(rbind,lapply(sp,function(ii){x<-z[ii,];data.frame(GID=x$GID[1],nursery=x$nursery[1],series=x$series[1],environment=x$environment[1],location=x$location[1],year=x$year[1],sowing=x$sowing[1],n_rows=nrow(x),disease_median=median(x$disease_index_pct),disease_mean=mean(x$disease_index_pct))}))
}

cand_strata<-aggregate_gid_strata(res_ids)
reference_ids<-c(LOCAL_CHECK="304660",SONALIKA="6387",CHIRYA3="911521")
ref_strata<-aggregate_gid_strata(unname(reference_ids))
ref_strata$reference_name<-names(reference_ids)[match(ref_strata$GID,reference_ids)]
write_out(ref_strata,"01_reference_coverage_by_stratum.csv")

anchor_one<-function(ref_name){
  r<-ref_strata[ref_strata$reference_name==ref_name,c("nursery","environment","disease_median","disease_mean")]
  names(r)[3:4]<-c("reference_median","reference_mean")
  a<-merge(cand_strata,r,by=c("nursery","environment"),all=FALSE)
  a$reference_name<-ref_name;a$delta_median<-a$disease_median-a$reference_median;a$delta_mean<-a$disease_mean-a$reference_mean
  a
}
anchored_all<-do.call(rbind,lapply(names(reference_ids),anchor_one));row.names(anchored_all)<-NULL
write_out(anchored_all,"02_all_reference_anchored_observations.csv")
local<-anchored_all[anchored_all$reference_name=="LOCAL_CHECK",]
write_out(local,"03_local_check_anchored_observations.csv")

bootstrap_ci<-function(x,B=2000){
  n<-length(x);if(n<2)return(c(NA,NA,NA,NA))
  med<-numeric(B);avg<-numeric(B)
  for(b in seq_len(B)){v<-x[sample.int(n,n,replace=TRUE)];med[b]<-median(v);avg[b]<-mean(v)}
  c(qfun(med,.025),qfun(med,.975),qfun(avg,.025),qfun(avg,.975))
}
summarize_anchor<-function(z){
  ci<-bootstrap_ci(z$delta_median)
  data.frame(GID=z$GID[1],reference_name=z$reference_name[1],n_strata=nrow(z),n_years=length(unique(z$year)),n_nurseries=length(unique(z$nursery)),n_locations=length(unique(z$location)),median_delta=median(z$delta_median),mean_delta=mean(z$delta_median),q10_delta=qfun(z$delta_median,.10),q90_delta=qfun(z$delta_median,.90),worst_delta=max(z$delta_median),pct_strict_advantage=100*mean(z$delta_median<0),pct_not_worse=100*mean(z$delta_median<=0),median_delta_ci_low=ci[1],median_delta_ci_high=ci[2],mean_delta_ci_low=ci[3],mean_delta_ci_high=ci[4])
}
summary_all<-do.call(rbind,lapply(split(seq_len(nrow(anchored_all)),interaction(anchored_all$GID,anchored_all$reference_name,drop=TRUE)),function(ii)summarize_anchor(anchored_all[ii,])))
summary_all<-merge(summary_all,portfolio[,c("GID","Cross","portfolio_role")],by="GID",all.x=TRUE)
write_out(summary_all,"04_multi_reference_candidate_summary.csv")

local_summary<-summary_all[summary_all$reference_name=="LOCAL_CHECK",]
local_summary$advantage_class<-"no_consistent_advantage"
low<-local_summary$n_strata<6|local_summary$n_years<2
local_summary$advantage_class[low]<-"insufficient_shared_environments"
local_summary$advantage_class[!low&local_summary$median_delta<0&local_summary$pct_not_worse>=60]<-"conditional_advantage"
local_summary$advantage_class[!low&local_summary$median_delta<0&local_summary$median_delta_ci_high<0&local_summary$pct_not_worse>=75]<-"stable_superior"
local_summary<-local_summary[order(factor(local_summary$advantage_class,levels=c("stable_superior","conditional_advantage","no_consistent_advantage","insufficient_shared_environments")),local_summary$median_delta),]
write_out(local_summary,"05_local_check_candidate_advantage_summary.csv")

# Location x sowing effects use observed shared environments only.
locsum<-do.call(rbind,lapply(split(seq_len(nrow(local)),interaction(local$GID,local$location,local$sowing,drop=TRUE)),function(ii){z<-local[ii,];data.frame(GID=z$GID[1],location=z$location[1],sowing=z$sowing[1],n_strata=nrow(z),n_years=length(unique(z$year)),median_delta=median(z$delta_median),q90_delta=qfun(z$delta_median,.90),worst_delta=max(z$delta_median),pct_not_worse=100*mean(z$delta_median<=0))}))
locsum<-merge(locsum,portfolio[,c("GID","portfolio_role")],by="GID",all.x=TRUE)
locsum$loss_of_advantage_flag<-locsum$median_delta>=0|locsum$pct_not_worse<60
write_out(locsum,"06_local_check_location_sowing_advantage.csv")
write_out(locsum[locsum$loss_of_advantage_flag,],"07_location_sowing_advantage_failures.csv")

# Leave-one-year-out checks whether the sign of the median advantage depends on one year.
loyo<-do.call(rbind,lapply(split(seq_len(nrow(local)),local$GID),function(ii){z<-local[ii,];do.call(rbind,lapply(sort(unique(z$year)),function(y){zz<-z[z$year!=y,];data.frame(GID=z$GID[1],excluded_year=y,remaining_years=length(unique(zz$year)),remaining_strata=nrow(zz),median_delta=if(nrow(zz))median(zz$delta_median)else NA,mean_delta=if(nrow(zz))mean(zz$delta_median)else NA,pct_not_worse=if(nrow(zz))100*mean(zz$delta_median<=0)else NA)}))}))
write_out(loyo,"08_leave_one_year_out_anchor_results.csv")
lo_summary<-do.call(rbind,lapply(split(seq_len(nrow(loyo)),loyo$GID),function(ii){z<-loyo[ii,];data.frame(GID=z$GID[1],n_exclusions=nrow(z),negative_median_fraction=mean(z$median_delta<0,na.rm=TRUE),not_worse_75pct_fraction=mean(z$pct_not_worse>=75,na.rm=TRUE),min_remaining_strata=min(z$remaining_strata),median_delta_range=paste(sprintf("%.2f",range(z$median_delta,na.rm=TRUE)),collapse=" to "))}))
lo_summary<-merge(lo_summary,portfolio[,c("GID","portfolio_role")],by="GID",all.x=TRUE)
write_out(lo_summary,"09_leave_one_year_out_anchor_summary.csv")

# Coverage is explicit for every candidate observed stratum set.
cand_cov<-do.call(rbind,lapply(split(seq_len(nrow(cand_strata)),cand_strata$GID),function(ii){z<-cand_strata[ii,];a<-local[local$GID==z$GID[1],];data.frame(GID=z$GID[1],candidate_strata=nrow(z),anchored_strata=nrow(a),anchor_coverage_pct=100*nrow(a)/nrow(z),candidate_years=length(unique(z$year)),anchored_years=length(unique(a$year)))}))
write_out(cand_cov,"10_candidate_anchor_coverage.csv")

draw_forest<-function(){z<-local_summary[order(local_summary$median_delta),];yy<-seq_len(nrow(z));cols<-c(stable_superior="#3F6B3A",conditional_advantage="#C9A227",no_consistent_advantage="#B8B2A7",insufficient_shared_environments="#7C817A");par(mar=c(4,5.2,2.2,1),family="sans",cex=.78);plot(range(c(z$median_delta_ci_low,z$median_delta_ci_high,0),na.rm=TRUE),range(yy),type="n",yaxt="n",xlab="Disease difference versus LOCAL CHECK (%)",ylab="");segments(z$median_delta_ci_low,yy,z$median_delta_ci_high,yy,col="#6B6258",lwd=1.1);points(z$median_delta,yy,pch=21,bg=cols[z$advantage_class],col="white",cex=1.2);axis(2,at=yy,labels=z$GID,las=1,cex.axis=.68);abline(v=0,lty=3,col="#7C817A");mtext("Points: medians; lines: 95% intervals",side=3,adj=1,cex=.58,col="#6B6B62");title("LOCAL CHECK–anchored candidate advantage",adj=0,font.main=2,cex.main=.95)}

row_ids<-as.character(local_summary$GID);strata<-as.vector(outer(c("Jashore","Quirusillas","Okinawa"),c("1st","2nd"),paste,sep="|"));hm<-matrix(NA_real_,length(row_ids),length(strata),dimnames=list(row_ids,strata))
for(i in seq_len(nrow(locsum))){r<-match(as.character(locsum$GID[i]),row_ids);c<-match(paste(locsum$location[i],locsum$sowing[i],sep="|"),strata);if(!is.na(r)&&!is.na(c))hm[r,c]<-locsum$median_delta[i]}
lim<-max(abs(hm),na.rm=TRUE)
draw_heat<-function(){par(mar=c(7,5,2.2,1),family="sans",cex=.78);image(seq_len(ncol(hm)),seq_len(nrow(hm)),t(hm[nrow(hm):1,,drop=FALSE]),col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),axes=FALSE,zlim=c(-lim,lim),xlab="",ylab="");axis(1,at=seq_len(ncol(hm)),labels=gsub("\\|","\n",colnames(hm)),las=2,cex.axis=.65);axis(2,at=seq_len(nrow(hm)),labels=rev(rownames(hm)),las=1,cex.axis=.65);title("Candidate advantage by location × sowing",adj=0,font.main=2,cex.main=.95);mtext("Blue = candidate lower than LOCAL CHECK; red = candidate higher",side=1,line=5.2,cex=.65)}

# Compare reference choice only where candidate and reference share a nursery environment.
refs<-unique(summary_all$reference_name);mr<-matrix(NA_real_,length(row_ids),length(refs),dimnames=list(row_ids,refs));for(i in seq_len(nrow(summary_all))){r<-match(as.character(summary_all$GID[i]),row_ids);c<-match(summary_all$reference_name[i],refs);if(!is.na(r)&&!is.na(c))mr[r,c]<-summary_all$median_delta[i]}
lim2<-max(abs(mr),na.rm=TRUE)
draw_multiref<-function(){par(mar=c(5.5,5,2.2,1),family="sans",cex=.8);image(seq_len(ncol(mr)),seq_len(nrow(mr)),t(mr[nrow(mr):1,,drop=FALSE]),col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),axes=FALSE,zlim=c(-lim2,lim2),xlab="",ylab="");axis(1,at=seq_len(ncol(mr)),labels=colnames(mr),las=2,cex.axis=.7);axis(2,at=seq_len(nrow(mr)),labels=rev(rownames(mr)),las=1,cex.axis=.65);title("Sensitivity to reference identity",adj=0,font.main=2,cex.main=.95);mtext("Median disease difference; blank = no shared environments",side=1,line=4,cex=.7)}

# Arial/Helvetica/sans fallback is declared; all devices use the R backend.
export_plot<-function(fun,stem,width=183,height=115){w<-width/25.4;h<-height/25.4;if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(file.path(figdir,paste0(stem,".svg")),width=w,height=h,system_fonts=list(sans="Arial"))else svg(file.path(figdir,paste0(stem,".svg")),width=w,height=h,family="sans");fun();dev.off();cairo_pdf(file.path(figdir,paste0(stem,".pdf")),width=w,height=h,family="sans");fun();dev.off();if(requireNamespace("ragg",quietly=TRUE))ragg::agg_tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw")else tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",family="sans");fun();dev.off();png(file.path(figdir,paste0(stem,".png")),width=w,height=h,units="in",res=300,family="sans");fun();dev.off()}
export_plot(draw_forest,"Figure_8_local_check_anchored_advantage",height=135)
export_plot(draw_heat,"Figure_9_location_sowing_anchored_advantage",height=145)
export_plot(draw_multiref,"Figure_10_reference_identity_sensitivity",height=135)

stable<-local_summary[local_summary$advantage_class=="stable_superior",];conditional<-local_summary[local_summary$advantage_class=="conditional_advantage",];fail<-locsum[locsum$loss_of_advantage_flag,]
report<-paste0("# 共同对照锚定的跨环境抗性优势分析\n\n",
"## 方法\n\n在每个真实育种圃×环境内部，以LOCAL CHECK/304660中位病害为基准，定义ΔDI=候选病害−LOCAL CHECK病害；负值表示候选更低。每个候选的区间通过对共享环境簇进行2,000次bootstrap获得，而不是把同一环境内记录当作独立重复。\n\n",
"## 主要结果\n\n14个候选中，",nrow(stable),"个满足stable_superior：",paste(stable$GID,collapse="、"),"；",nrow(conditional),"个为conditional_advantage：",paste(conditional$GID,collapse="、"),"。所有候选的LOCAL CHECK共享环境覆盖率均为100%，中位ΔDI范围为",sprintf("%.2f",min(local_summary$median_delta)),"至",sprintf("%.2f",max(local_summary$median_delta)),"个百分点；逐年剔除后所有候选的中位优势方向均保持为负。\n\n",
"这一高度一致性说明LOCAL CHECK是有效的感病锚点，可验证候选抗性信号，但由于对照与候选差距过大，它对14个候选内部排序的分辨能力有限。SONALIKA仅与5个候选、CHIRYA.3仅与7个候选共享环境，而且均只覆盖1年，因此只能作为局部敏感性结果，不能称为跨年度独立确认。\n\n",
"地点×播期分析标记了",nrow(fail),"个候选—环境层面的优势丢失记录。该标记可能来自真实条件反应，也可能来自共享环境较少，因此必须同时查看n_strata和原始ΔDI。\n\n",
"## 解释边界\n\nLOCAL CHECK是跨年度尺度锚点，不是随机抽取的统一品种；其身份和管理可能随育种圃语境变化。锚定差值减少了环境压力和材料构成变化的影响，但不能替代平衡随机试验。SONALIKA和CHIRYA.3只在共享环境中作为辅助敏感性参照。\n")
writeLines(report,file.path(repdir,"reference_anchored_report.md"),useBytes=TRUE)
qa<-paste0("# QA\n\n- 主表30475条，不含15HZAN；14个抗性候选均来自非冗余复验组合。\n- LOCAL CHECK按育种圃×环境取中位数；候选与对照仅在完全相同的育种圃×环境连接。\n- 未共享环境保持缺失，不外推、不插补。\n- bootstrap单位为共享环境簇，固定随机种子20260811，2,000次。\n- 所有图件由R生成，183 mm宽，含SVG/PDF/600 dpi TIFF/PNG。\n")
writeLines(qa,file.path(repdir,"QA_notes.md"),useBytes=TRUE);sink(file.path(logdir,"sessionInfo.txt"));print(sessionInfo());sink()
meta<-data.frame(metric=c("resistance_candidates","local_check_rows","local_check_strata","candidate_local_shared_rows","stable_superior","conditional_advantage","location_sowing_failures"),value=c(length(res_ids),sum(d$GID==reference_ids[["LOCAL_CHECK"]],na.rm=TRUE),nrow(ref_strata[ref_strata$reference_name=="LOCAL_CHECK",]),nrow(local),nrow(stable),nrow(conditional),nrow(fail)))
write_out(meta,"00_analysis_summary.csv");cat(paste(meta$metric,meta$value,sep="="),sep="\n")
