options(stringsAsFactors=FALSE,warn=1)
set.seed(20260811)
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);base<-Sys.getenv("CIMMYT_ANALYSIS_ROOT",unset=root)
longfile<-file.path(base,"history","standardized","historical_long_primary_observed.csv")
portfile<-file.path(base,"history","validation_portfolio_analysis","tables","01_nonredundant_validation_portfolio.csv")
out<-file.path(base,"history","temporal_validation_analysis");tabdir<-file.path(out,"tables");figdir<-file.path(out,"figures");repdir<-file.path(out,"reports");logdir<-file.path(out,"logs")
for(p in c(tabdir,figdir,repdir,logdir))dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_out<-function(x,n)write.csv(x,file.path(tabdir,n),row.names=FALSE,na="",fileEncoding="UTF-8")
qfun<-function(x,p)unname(quantile(x,p,type=7,na.rm=TRUE));safe_cor<-function(x,y)if(length(x)>=3&&sd(x)>0&&sd(y)>0)unname(cor(x,y,method="spearman"))else NA_real_
d<-read.csv(longfile,check.names=FALSE,na.strings=c("","NA"),encoding="UTF-8-BOM");portfolio<-read.csv(portfile,check.names=FALSE)
d$GID<-trimws(as.character(d$GID));d$GID[d$GID==""]<-NA_character_;d$year<-as.integer(d$year);d$disease_index_pct<-as.numeric(d$disease_index_pct)
stopifnot(nrow(d)==30475,!any(d$nursery=="15HZAN"));d<-d[!is.na(d$GID),]

# Scores are created within each observed nursery x environment and never use future years.
g_ne<-interaction(d$nursery,d$environment,drop=TRUE)
d$within_environment_percentile<-ave(d$disease_index_pct,g_ne,FUN=function(x)rank(x,ties.method="average")/length(x))
check<-d[d$GID=="304660",]
check_med<-do.call(rbind,lapply(split(seq_len(nrow(check)),interaction(check$nursery,check$environment,drop=TRUE)),function(ii){z<-check[ii,];data.frame(nursery=z$nursery[1],environment=z$environment[1],check_median=median(z$disease_index_pct))}))
d<-merge(d,check_med,by=c("nursery","environment"),all.x=TRUE,sort=FALSE)
d$delta_local_check<-d$disease_index_pct-d$check_median

summarize_rows<-function(z,prefix=""){
  data.frame(n_obs=nrow(z),n_years=length(unique(z$year)),n_nurseries=length(unique(z$nursery)),n_locations=length(unique(z$location)),n_environments=length(unique(paste(z$nursery,z$environment))),median_percentile=median(z$within_environment_percentile),q75_percentile=qfun(z$within_environment_percentile,.75),q90_disease=qfun(z$disease_index_pct,.90),worst_disease=max(z$disease_index_pct),median_delta=if(all(is.na(z$delta_local_check)))NA_real_ else median(z$delta_local_check,na.rm=TRUE),pct_not_worse_check=if(all(is.na(z$delta_local_check)))NA_real_ else 100*mean(z$delta_local_check<=0,na.rm=TRUE))
}

gy<-do.call(rbind,lapply(split(seq_len(nrow(d)),interaction(d$GID,d$year,drop=TRUE)),function(ii){z<-d[ii,];cbind(data.frame(GID=z$GID[1],year=z$year[1],Cross=names(sort(table(z$Cross),decreasing=TRUE))[1]),summarize_rows(z))}))
row.names(gy)<-NULL;write_out(gy,"01_gid_year_observed_metrics.csv")

discovery_rule<-function(z){z$n_obs>=3&z$n_locations>=2&z$median_percentile<=.40&z$q75_percentile<=.50&z$q90_disease<=20&(is.na(z$median_delta)|z$median_delta<0)}
test_rule<-function(z){z$n_obs>=3&z$median_percentile<=.50&z$q75_percentile<=.65&z$q90_disease<=30&(is.na(z$median_delta)|z$median_delta<0)}

transition_rows<-list();transition_summary<-list();counter<-1
for(cut in 2018:2023){
  test_year<-cut+1;train_ids<-unique(d$GID[d$year<=cut]);test_ids<-unique(d$GID[d$year==test_year]);common<-setdiff(intersect(train_ids,test_ids),"304660")
  rows<-list();ri<-1
  for(id in common){
    tr<-d[d$GID==id&d$year<=cut,];te<-d[d$GID==id&d$year==test_year,]
    a<-summarize_rows(tr);b<-summarize_rows(te)
    rows[[ri]]<-data.frame(cutoff_year=cut,test_year=test_year,GID=id,train_n_obs=a$n_obs,train_n_years=a$n_years,train_n_locations=a$n_locations,train_n_environments=a$n_environments,train_median_percentile=a$median_percentile,train_q75_percentile=a$q75_percentile,train_q90_disease=a$q90_disease,train_median_delta=a$median_delta,test_n_obs=b$n_obs,test_n_locations=b$n_locations,test_median_percentile=b$median_percentile,test_q75_percentile=b$q75_percentile,test_q90_disease=b$q90_disease,test_median_delta=b$median_delta)
    ri<-ri+1
  }
  trn<-if(length(rows))do.call(rbind,rows)else data.frame()
  if(nrow(trn)){
    trn$eligible_train<-trn$train_n_obs>=3&trn$train_n_locations>=2&trn$test_n_obs>=3
    trn$selected_train<-FALSE;trn$selected_train[trn$eligible_train]<-discovery_rule(data.frame(n_obs=trn$train_n_obs[trn$eligible_train],n_locations=trn$train_n_locations[trn$eligible_train],median_percentile=trn$train_median_percentile[trn$eligible_train],q75_percentile=trn$train_q75_percentile[trn$eligible_train],q90_disease=trn$train_q90_disease[trn$eligible_train],median_delta=trn$train_median_delta[trn$eligible_train]))
    trn$passed_test<-test_rule(data.frame(n_obs=trn$test_n_obs,n_locations=trn$test_n_locations,median_percentile=trn$test_median_percentile,q75_percentile=trn$test_q75_percentile,q90_disease=trn$test_q90_disease,median_delta=trn$test_median_delta))
    eligible<-trn[trn$eligible_train,]
    if(nrow(eligible)){
      eligible$train_rank_fraction<-rank(eligible$train_median_percentile,ties.method="average")/nrow(eligible)
      eligible$test_rank_fraction<-rank(eligible$test_median_percentile,ties.method="average")/nrow(eligible)
      top_train<-eligible$train_rank_fraction<=.25;top_test<-eligible$test_rank_fraction<=.25
      precision<-if(sum(top_train))mean(top_test[top_train])else NA_real_
      recall<-if(sum(top_test))mean(top_train[top_test])else NA_real_
      rho<-safe_cor(eligible$train_median_percentile,eligible$test_median_percentile)
      selected_n<-sum(eligible$selected_train);selected_precision<-if(selected_n)mean(eligible$passed_test[eligible$selected_train])else NA_real_
      transition_summary[[counter]]<-data.frame(cutoff_year=cut,test_year=test_year,shared_gid=length(common),eligible_gid=nrow(eligible),selected_gid=selected_n,selected_test_pass=if(selected_n)sum(eligible$passed_test[eligible$selected_train])else 0,selected_precision=selected_precision,top_quartile_precision=precision,top_quartile_recall=recall,spearman_train_test=rho)
      trn<-merge(trn,eligible[,c("GID","train_rank_fraction","test_rank_fraction")],by="GID",all.x=TRUE)
    }else transition_summary[[counter]]<-data.frame(cutoff_year=cut,test_year=test_year,shared_gid=length(common),eligible_gid=0,selected_gid=0,selected_test_pass=0,selected_precision=NA,top_quartile_precision=NA,top_quartile_recall=NA,spearman_train_test=NA)
  }else transition_summary[[counter]]<-data.frame(cutoff_year=cut,test_year=test_year,shared_gid=0,eligible_gid=0,selected_gid=0,selected_test_pass=0,selected_precision=NA,top_quartile_precision=NA,top_quartile_recall=NA,spearman_train_test=NA)
  transition_rows[[counter]]<-trn;counter<-counter+1
}
pred<-do.call(rbind,transition_rows);tsum<-do.call(rbind,transition_summary)
write_out(pred,"02_all_temporal_transition_predictions.csv");write_out(tsum,"03_temporal_transition_summary.csv")

# Aggregate prediction evidence across all material-transition events.
valid_sel<-pred[pred$eligible_train&pred$selected_train,]
overall<-data.frame(metric=c("transitions_with_eligible_gid","eligible_gid_transition_events","selected_transition_events","passed_selected_events","pooled_selected_precision","median_transition_precision","median_top_quartile_precision","median_spearman"),value=c(sum(tsum$eligible_gid>0),sum(tsum$eligible_gid),nrow(valid_sel),sum(valid_sel$passed_test),if(nrow(valid_sel))mean(valid_sel$passed_test)else NA,median(tsum$selected_precision,na.rm=TRUE),median(tsum$top_quartile_precision,na.rm=TRUE),median(tsum$spearman_train_test,na.rm=TRUE)))
write_out(overall,"04_overall_temporal_validation_metrics.csv")

# Retrospective chronology for the fixed 14-material portfolio. This is not independent selection.
port_ids<-as.character(portfolio$GID[portfolio$selected_for_resistance_validation]);port_gy<-gy[gy$GID%in%port_ids,]
port_gy$year_test_pass<-test_rule(port_gy)
first_year<-tapply(port_gy$year,port_gy$GID,min)
port_gy$chronology_role<-ifelse(port_gy$year==first_year[port_gy$GID],"first_observed_discovery","later_observed_test")
write_out(port_gy,"05_portfolio_chronological_year_outcomes.csv")
port_summary<-do.call(rbind,lapply(split(seq_len(nrow(port_gy)),port_gy$GID),function(ii){z<-port_gy[ii,];later<-z[z$chronology_role=="later_observed_test",];data.frame(GID=z$GID[1],first_year=min(z$year),last_year=max(z$year),n_years=nrow(z),first_year_pass=z$year_test_pass[z$year==min(z$year)][1],later_years=nrow(later),later_pass_fraction=if(nrow(later))mean(later$year_test_pass)else NA,later_median_percentile=if(nrow(later))median(later$median_percentile)else NA,later_q90_disease_max=if(nrow(later))max(later$q90_disease)else NA)}))
port_summary<-merge(port_summary,portfolio[,c("GID","portfolio_role")],by="GID",all.x=TRUE);write_out(port_summary,"06_portfolio_temporal_validation_summary.csv")

# Three pre-specified discovery rules illustrate selection-performance trade-offs.
scenarios<-data.frame(scenario=c("conservative","reference","liberal"),median_pct=c(.30,.40,.50),q75_pct=c(.40,.50,.65),q90_disease=c(10,20,30))
scenario_eval<-do.call(rbind,lapply(seq_len(nrow(scenarios)),function(i){s<-scenarios[i,];z<-pred[pred$eligible_train,];sel<-z$train_median_percentile<=s$median_pct&z$train_q75_percentile<=s$q75_pct&z$train_q90_disease<=s$q90_disease&(is.na(z$train_median_delta)|z$train_median_delta<0);data.frame(scenario=s$scenario,selected_events=sum(sel),passed_events=sum(z$passed_test[sel]),pooled_precision=if(sum(sel))mean(z$passed_test[sel])else NA,distinct_gid=length(unique(z$GID[sel])))}))
write_out(scenario_eval,"07_temporal_rule_sensitivity.csv")

draw_transitions<-function(){layout(matrix(c(1,2),1,2));par(mar=c(6.8,4.5,2.2,.8),mgp=c(2.5,.7,0),family="sans",cex=.8);labs<-paste0(tsum$cutoff_year,"→",tsum$test_year);barplot(tsum$selected_precision,names.arg=labs,ylim=c(0,1),col="#587A45",border=NA,las=2,ylab="Conditional next-year pass proportion");abline(h=.5,lty=3,col="#666666");title("a  Conditional next-year performance",adj=0,font.main=2,cex.main=.95);plot(seq_len(nrow(tsum)),tsum$spearman_train_test,type="b",pch=21,bg="#A66A3F",col="#A66A3F",ylim=c(-1,1),xaxt="n",xlab="",ylab="Train–test Spearman correlation");axis(1,at=seq_len(nrow(tsum)),labels=labs,las=2,cex.axis=.65);abline(h=0,lty=3,col="#666666");text(seq_len(nrow(tsum)),tsum$spearman_train_test,labels=paste0("n=",tsum$eligible_gid),pos=3,cex=.55);title("b  Rank transportability",adj=0,font.main=2,cex.main=.95)}

valid<-pred[pred$eligible_train,];draw_scatter<-function(){cuts<-sort(unique(valid$cutoff_year));layout(matrix(seq_len(6),2,3,byrow=TRUE));par(mar=c(3,3,1.8,.5),mgp=c(1.8,.55,0),tcl=-.2,family="sans",cex=.7);for(cut in 2018:2023){z<-valid[valid$cutoff_year==cut,];if(nrow(z)){plot(z$train_median_percentile,z$test_median_percentile,pch=ifelse(z$selected_train,21,1),bg=ifelse(z$selected_train,"#228833","white"),col="#555555",xlim=c(0,1),ylim=c(0,1),xlab="Train median percentile",ylab="Next-year median percentile");abline(0,1,lty=3,col="#999999");abline(v=.4,h=.5,lty=3,col="#777777");title(paste0(cut,"→",cut+1," (n=",nrow(z),")"),adj=0,font.main=2,cex.main=.9)}else{plot.new();title(paste0(cut,"→",cut+1," (no eligible GID)"),adj=0,font.main=2,cex.main=.9)}}}

years<-2018:2024;hm<-matrix(NA_real_,length(port_ids),length(years),dimnames=list(port_ids,years));for(i in seq_len(nrow(port_gy))){r<-match(port_gy$GID[i],port_ids);c<-match(port_gy$year[i],years);hm[r,c]<-ifelse(port_gy$year_test_pass[i],1,0)}
draw_portfolio<-function(){par(mar=c(4,5,2.2,1),family="sans",cex=.8);image(seq_along(years),seq_along(port_ids),t(hm[length(port_ids):1,,drop=FALSE]),col=c("#CC6677","#228833"),axes=FALSE,zlim=c(0,1),xlab="Year",ylab="");axis(1,at=seq_along(years),labels=years);axis(2,at=seq_along(port_ids),labels=rev(port_ids),las=1,cex.axis=.65);title("Retrospective yearly outcomes of the fixed portfolio",adj=0,font.main=2,cex.main=.95);legend("bottomright",legend=c("test rule failed","test rule passed","not observed"),fill=c("#CC6677","#228833","white"),bty="n",cex=.65)}

# Arial/Helvetica/sans fallback is declared; all devices use the R backend.
export_plot<-function(fun,stem,width=183,height=110){w<-width/25.4;h<-height/25.4;if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(file.path(figdir,paste0(stem,".svg")),width=w,height=h,system_fonts=list(sans="Arial"))else svg(file.path(figdir,paste0(stem,".svg")),width=w,height=h,family="sans");fun();dev.off();cairo_pdf(file.path(figdir,paste0(stem,".pdf")),width=w,height=h,family="sans");fun();dev.off();if(requireNamespace("ragg",quietly=TRUE))ragg::agg_tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw")else tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",family="sans");fun();dev.off();png(file.path(figdir,paste0(stem,".png")),width=w,height=h,units="in",res=300,family="sans");fun();dev.off()}
export_plot(draw_transitions,"Figure_11_temporal_selection_and_rank_transport")
export_plot(draw_scatter,"Figure_12_train_next_year_scatter",height=125)
export_plot(draw_portfolio,"Figure_13_portfolio_yearly_outcomes",height=125)

best<-tsum[which.max(tsum$selected_precision),];worst<-tsum[which.min(tsum$selected_precision),]
report<-paste0("# 时间外伪前瞻验证\n\n",
"## 设计\n\n每个年度转移仅使用截至年份t的数据计算训练表现，并用t+1年度评价。主分析包含所有训练期和测试期均有真实GID连接的材料；训练规则要求至少3个观测、2个地点、中位百分位≤0.40、75分位≤0.50、病害90分位≤20，并在可锚定时保持LOCAL CHECK优势。测试通过要求下一年中位百分位≤0.50、75分位≤0.65、病害90分位≤30且对照优势不反转。\n\n",
"## 结果\n\n共有",sum(tsum$eligible_gid>0),"个年度转移具有可评价材料，合计",sum(tsum$eligible_gid),"个合格GID—转移事件。训练规则选出",nrow(valid_sel),"个事件，其中",sum(valid_sel$passed_test),"个在下一年度通过，合并精确率为",sprintf("%.1f",100*mean(valid_sel$passed_test)),"%。年度间精确率范围为",sprintf("%.1f",100*min(tsum$selected_precision,na.rm=TRUE)),"%至",sprintf("%.1f",100*max(tsum$selected_precision,na.rm=TRUE)),"%。\n\n",
"固定14材料组合的逐年结果属于回溯性轨迹，不是独立选择验证；其结果见`06_portfolio_temporal_validation_summary.csv`。环境与材料组成在年度间同时改变，因此排名相关性衡量可迁移性，而不是遗传决定性。\n\n",
"## 限制\n\n只有在相邻年度重复出现的GID才能进入测试，样本是非随机存留子集；部分晚期转移可能只有少量对照或长期材料。零膨胀会压缩秩差异。该分析减少时间泄漏，但仍是历史伪前瞻验证，不能替代真正前瞻性试验。\n")
writeLines(report,file.path(repdir,"temporal_validation_report.md"),useBytes=TRUE)
qa<-paste0("# QA\n\n- 输入30475条主分析记录，不含15HZAN。\n- 所有训练统计仅使用year≤cutoff；测试统计仅使用下一年度。\n- 主分析覆盖全部可连接GID，不以最终14候选限定样本。\n- 不构造缺失材料×年份单元；无重复GID的转移保持不可评价。\n- 固定组合轨迹明确标记为回溯性，不称为独立验证。\n- 图件由R生成，183 mm宽，SVG/PDF/600 dpi TIFF/PNG。\n")
writeLines(qa,file.path(repdir,"QA_notes.md"),useBytes=TRUE);sink(file.path(logdir,"sessionInfo.txt"));print(sessionInfo());sink()
meta<-data.frame(metric=c("transitions","transitions_evaluable","eligible_transition_events","selected_transition_events","passed_selected_events","pooled_precision","portfolio_gid"),value=c(6,sum(tsum$eligible_gid>0),sum(tsum$eligible_gid),nrow(valid_sel),sum(valid_sel$passed_test),if(nrow(valid_sel))mean(valid_sel$passed_test)else NA,length(port_ids)))
write_out(meta,"00_analysis_summary.csv");cat(paste(meta$metric,meta$value,sep="="),sep="\n")
