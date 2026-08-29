options(stringsAsFactors=FALSE, scipen=999)
set.seed(20260811)

root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
base <- file.path(root,"CIMMYT数据集")
tab <- file.path(base,"结果表"); fig <- file.path(base,"图表"); repd <- file.path(base,"报告"); logd <- file.path(base,"运行记录")
for(d in c(tab,fig,repd,logd)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
input_obs <- file.path(root,"外部验证数据","CIMMYT_2024","处理后数据","cimmyt_wheat_blast_long_observed.csv")
input_all <- file.path(root,"外部验证数据","CIMMYT_2024","处理后数据","cimmyt_wheat_blast_long_all.csv")
phase1 <- file.path(tab,"03_材料跨环境表现与候选筛选.csv")
stopifnot(file.exists(input_obs),file.exists(input_all),file.exists(phase1))

obs <- read.csv(input_obs,check.names=FALSE,fileEncoding="UTF-8-BOM")
all_dat <- read.csv(input_all,check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings=c(""))
mat1 <- read.csv(phase1,check.names=FALSE,fileEncoding="UTF-8-BOM")
for(nm in c("obs","all_dat")) {
  z <- get(nm); names(z)[names(z)=="Ent"] <- "Entry"; z$Entry <- as.character(z$Entry); z$GID <- as.character(z$GID)
  z$disease_index_pct <- suppressWarnings(as.numeric(z$disease_index_pct)); z$material_id <- paste(z$nursery,z$Entry,sep="::")
  z$env_id <- paste(z$location,z$year,z$sowing,sep="_"); assign(nm,z)
}
nurseries <- sort(unique(obs$nursery))
qv <- function(x,p) as.numeric(quantile(x,p,names=FALSE,type=7,na.rm=TRUE))
wcsv <- function(x,name) write.csv(x,file.path(tab,name),row.names=FALSE,na="",fileEncoding="UTF-8")
fmt <- function(x,d=2) ifelse(is.na(x),"NA",format(round(x,d),nsmall=d,trim=TRUE))

# A. Environment discrimination: four bounded components and bootstrap uncertainty.
entropy5 <- function(x) {
  g <- ifelse(x==0,"zero",ifelse(x<=5,"gt0_5",ifelse(x<=10,"gt5_10",ifelse(x<20,"gt10_20","ge20"))))
  p <- table(factor(g,levels=c("zero","gt0_5","gt5_10","gt10_20","ge20")))/length(x); p <- p[p>0]
  -sum(p*log(p))/log(5)
}
disc_metrics <- function(x) {
  nonzero <- mean(x>0); iq <- IQR(x); spread <- qv(x,.9)-qv(x,.1); ent <- entropy5(x)
  comps <- c(nonzero_component=pmin(nonzero/.50,1),iqr_component=pmin(iq/10,1),spread_component=pmin(spread/20,1),entropy_component=ent)
  c(n=length(x),mean=mean(x),median=median(x),zero_prop=mean(x==0),sd=sd(x),mad=mad(x),iqr=iq,q10=qv(x,.1),q90=qv(x,.9),
    ge20_prop=mean(x>=20),normalized_entropy=ent,comps,discrimination_score=100*mean(comps))
}
env_groups <- split(obs,interaction(obs$nursery,obs$environment,drop=TRUE))
env_rows <- lapply(env_groups,function(d) {
  m <- disc_metrics(d$disease_index_pct)
  data.frame(nursery=d$nursery[1],environment=d$environment[1],env_id=d$env_id[1],location=d$location[1],year=d$year[1],sowing=d$sowing[1],t(m),check.names=FALSE)
})
env_disc <- do.call(rbind,env_rows)
# Bootstrap samples material rows within each nursery-environment; descriptive uncertainty only.
ci_rows <- lapply(env_groups,function(d) {
  x <- d$disease_index_pct; bs <- replicate(1000,disc_metrics(sample(x,length(x),replace=TRUE))["discrimination_score"])
  data.frame(nursery=d$nursery[1],environment=d$environment[1],score_ci_low=qv(bs,.025),score_ci_high=qv(bs,.975))
})
env_disc <- merge(env_disc,do.call(rbind,ci_rows),by=c("nursery","environment"),all.x=TRUE)
env_disc$score_class <- cut(env_disc$discrimination_score,c(-Inf,15,35,60,Inf),labels=c("minimal","low","moderate","high"),right=FALSE)
env_disc$rule_A <- env_disc$zero_prop<.80 & env_disc$iqr>0
env_disc$rule_B <- env_disc$zero_prop<.90 & env_disc$iqr>=5
env_disc$rule_C <- env_disc$discrimination_score>=35
env_disc$rule_count <- rowSums(env_disc[c("rule_A","rule_B","rule_C")])
env_disc$consensus_discriminating <- env_disc$rule_count>=2
env_disc <- env_disc[order(env_disc$nursery,-env_disc$discrimination_score),]
wcsv(env_disc,"15_环境区分力综合评价.csv")

disc_key <- paste(env_disc$nursery[env_disc$consensus_discriminating],env_disc$environment[env_disc$consensus_discriminating],sep="::")
obs$disc_env <- paste(obs$nursery,obs$environment,sep="::") %in% disc_key
obs$relative_resistance <- NA_real_
for(idx in split(seq_len(nrow(obs)),interaction(obs$nursery,obs$environment,drop=TRUE))) {
  x <- obs$disease_index_pct[idx]; obs$relative_resistance[idx] <- if(length(x)==1) 1 else 1-(rank(x,ties.method="average")-1)/(length(x)-1)
}

# B. Candidate robustness across six transparent selection rules.
rob_rows <- lapply(split(obs,obs$material_id),function(d) {
  di <- d[d$disc_env,]; x <- d$disease_index_pct; xi <- di$disease_index_pct; ri <- di$relative_resistance
  p1 <- mat1[mat1$material_id==d$material_id[1],]
  ndisc_nur <- sum(env_disc$nursery==d$nursery[1]&env_disc$consensus_discriminating)
  eligible <- ndisc_nur>=2 & length(xi)>=2 & p1$complete_rate[1]>=.8
  s1 <- eligible && mean(xi)<=5 && max(xi)<=10 && mean(xi<=5)>=.8
  s2 <- eligible && mean(xi)<=10 && max(xi)<=20 && mean(xi<=10)>=.8
  s3 <- eligible && mean(ri)>=.80 && min(ri)>=.50
  s4 <- eligible && mean(x)<=5 && max(x)<=10
  # S5 assigned below from within-nursery composite rank.
  s6 <- eligible && mean(xi<=5)>=.8 && !any(xi>=20)
  data.frame(nursery=d$nursery[1],Entry=d$Entry[1],GID=d$GID[1],Cross=d$Cross[1],Sel_Hist=d$Sel_Hist[1],material_id=d$material_id[1],
    eligible=eligible,n_discriminating_env=length(xi),disc_mean=if(length(xi))mean(xi)else NA,disc_max=if(length(xi))max(xi)else NA,
    disc_low5_prop=if(length(xi))mean(xi<=5)else NA,disc_mean_relative=if(length(ri))mean(ri)else NA,all_env_mean=mean(x),all_env_max=max(x),
    S1_strict_raw=s1,S2_relaxed_raw=s2,S3_relative_rank=s3,S4_all_env_worst_case=s4,S6_no_severe_failure=s6,
    phase1_composite=p1$composite_score[1])
})
rob <- do.call(rbind,rob_rows)
rob$S5_top20_composite <- FALSE
for(nur in nurseries) {
  ix <- which(rob$nursery==nur&rob$eligible); if(length(ix)) {
    cutoff <- qv(rob$phase1_composite[ix],.80); rob$S5_top20_composite[ix] <- rob$phase1_composite[ix]>=cutoff
  }
}
sc <- c("S1_strict_raw","S2_relaxed_raw","S3_relative_rank","S4_all_env_worst_case","S5_top20_composite","S6_no_severe_failure")
rob$selection_count <- rowSums(rob[sc]); rob$selection_frequency <- rob$selection_count/length(sc)
rob$robustness_class <- ifelse(!rob$eligible,"not_screenable",ifelse(rob$selection_count>=4,"robust",ifelse(rob$selection_count>=2,"provisional","weak")))
rob <- rob[order(rob$nursery,-rob$selection_count,rob$disc_mean),]
wcsv(rob,"16_候选材料六规则稳健性矩阵.csv")
rob_sum <- do.call(rbind,lapply(split(rob,rob$nursery),function(d) {
  data.frame(nursery=d$nursery[1],n_materials=nrow(d),eligible=sum(d$eligible),robust=sum(d$robustness_class=="robust"),
    provisional=sum(d$robustness_class=="provisional"),weak=sum(d$robustness_class=="weak"),not_screenable=sum(d$robustness_class=="not_screenable"))
}))
wcsv(rob_sum,"17_候选稳健性数量汇总.csv")

# C. Pairwise rank reversal in consensus-discriminating environments.
rank_pairs <- list(); rp <- 0L
for(nur in nurseries) {
  d <- obs[obs$nursery==nur&obs$disc_env,]; es <- sort(unique(d$env_id)); if(length(es)<2) next
  for(i in 1:(length(es)-1)) for(j in (i+1):length(es)) {
    a<-d[d$env_id==es[i],c("material_id","disease_index_pct","relative_resistance")]; b<-d[d$env_id==es[j],c("material_id","disease_index_pct","relative_resistance")]
    names(a)[2:3]<-c("x","rx"); names(b)[2:3]<-c("y","ry"); m<-merge(a,b,by="material_id")
    # Raw-scale quintile thresholds retain all ties; sets can exceed exactly 20% in zero-inflated data.
    topx<-m$x<=qv(m$x,.20); topy<-m$y<=qv(m$y,.20); botx<-m$x>=qv(m$x,.80); boty<-m$y>=qv(m$y,.80)
    rp<-rp+1L; rank_pairs[[rp]]<-data.frame(nursery=nur,environment1=es[i],environment2=es[j],n_pairs=nrow(m),
      spearman_rho=if(sd(m$x)>0&&sd(m$y)>0)cor(m$x,m$y,method="spearman")else NA,
      kendall_tau=if(sd(m$x)>0&&sd(m$y)>0)cor(m$x,m$y,method="kendall")else NA,
      mean_abs_percentile_change=mean(abs(m$rx-m$ry)),max_abs_percentile_change=max(abs(m$rx-m$ry)),
      top20_overlap_n=sum(topx&topy),top20_jaccard=sum(topx&topy)/max(1,sum(topx|topy)),
      top_to_bottom_n=sum((topx&boty)|(topy&botx)),top_to_bottom_prop=mean((topx&boty)|(topy&botx)))
  }
}
rank_pair_df <- if(length(rank_pairs))do.call(rbind,rank_pairs)else data.frame()
wcsv(rank_pair_df,"18_环境间排名翻转汇总.csv")

rank_mat <- lapply(split(obs[obs$disc_env,],obs$material_id[obs$disc_env]),function(d) {
  r<-d$relative_resistance
  cls<-if(length(r)<2)"insufficient" else if(mean(r)>=.75&&min(r)>=.50)"consistently_low_disease_rank" else if(mean(r)<=.25&&max(r)<=.50)"consistently_high_disease_rank" else if(diff(range(r))>=.50)"environment_specific" else "intermediate"
  data.frame(nursery=d$nursery[1],Entry=d$Entry[1],GID=d$GID[1],material_id=d$material_id[1],n_env=length(r),mean_relative=mean(r),min_relative=min(r),max_relative=max(r),rank_range=diff(range(r)),rank_sd=if(length(r)>1)sd(r)else NA,top20_prop=mean(r>=.8),bottom20_prop=mean(r<=.2),rank_pattern=cls)
})
rank_mat_df <- if(length(rank_mat))do.call(rbind,rank_mat)else data.frame()
wcsv(rank_mat_df,"19_材料排名稳定性与环境特异性.csv")

# D. Material-level sowing responses, preserving nursery-location-year pairing.
sow_pairs <- list(); spn <- 0L
for(k in split(obs,interaction(obs$nursery,obs$location,obs$year,drop=TRUE))) {
  a<-k[k$sowing=="1st",c("nursery","location","year","material_id","Entry","GID","disease_index_pct")]
  b<-k[k$sowing=="2nd",c("material_id","disease_index_pct")]; if(!nrow(a)||!nrow(b))next
  names(a)[7]<-"first"; names(b)[2]<-"second"; m<-merge(a,b,by="material_id"); m$delta<-m$second-m$first; m$abs_delta<-abs(m$delta)
  m$response_class<-ifelse(m$first<=5&m$second<=5,"both_low",ifelse(m$first<=5&m$second>10,"second_sowing_increase",ifelse(m$first>10&m$second<=5,"second_sowing_decrease",ifelse(m$first>10&m$second>10,"both_high","moderate_or_mixed"))))
  spn<-spn+1L;sow_pairs[[spn]]<-m
}
sow_pair_df<-do.call(rbind,sow_pairs); wcsv(sow_pair_df,"20_材料地点层面播期响应.csv")
sow_mat <- lapply(split(sow_pair_df,sow_pair_df$material_id),function(d) {
  inc<-any(d$response_class=="second_sowing_increase"); dec<-any(d$response_class=="second_sowing_decrease")
  cls<-if(all(d$response_class=="both_low"))"stable_low_both_sowings" else if(inc&&dec)"bidirectional_sensitive" else if(inc)"second_sowing_sensitive_increase" else if(dec)"second_sowing_sensitive_decrease" else if(mean(d$abs_delta<=5)>=.8)"generally_stable" else "variable_moderate"
  data.frame(nursery=d$nursery[1],Entry=d$Entry[1],GID=d$GID[1],material_id=d$material_id[1],n_location_year_pairs=nrow(d),mean_delta=mean(d$delta),median_delta=median(d$delta),mean_abs_delta=mean(d$abs_delta),max_increase=max(d$delta),max_decrease=min(d$delta),prop_abs_delta_le5=mean(d$abs_delta<=5),sowing_response_pattern=cls)
})
sow_mat_df<-do.call(rbind,sow_mat)
sow_mat_df<-merge(sow_mat_df,rob[c("material_id","robustness_class","selection_count")],by="material_id",all.x=TRUE)
sow_mat_df<-sow_mat_df[order(sow_mat_df$nursery,sow_mat_df$sowing_response_pattern,-sow_mat_df$selection_count),]
wcsv(sow_mat_df,"21_材料播期响应综合分类.csv")
sow_cross <- as.data.frame.matrix(table(sow_mat_df$nursery,sow_mat_df$sowing_response_pattern)); sow_cross$nursery<-rownames(sow_cross); rownames(sow_cross)<-NULL; sow_cross<-sow_cross[c("nursery",setdiff(names(sow_cross),"nursery"))]
wcsv(sow_cross,"22_播期响应类型数量汇总.csv")

# Figures
class_col<-c(minimal="#bdbdbd",low="#fee08b",moderate="#fc8d59",high="#d73027")
png(file.path(fig,"图6_环境区分力评分.png"),width=1700,height=1050,res=170)
par(mar=c(9,5,3,1)); z<-env_disc[order(env_disc$discrimination_score),]; labs<-paste0(z$nursery,"\n",substr(z$location,1,3),substr(z$year,3,4),"-S",substr(z$sowing,1,1))
bp<-barplot(z$discrimination_score,names.arg=labs,las=2,col=class_col[as.character(z$score_class)],border=NA,ylim=c(0,100),ylab="Discrimination score (0-100)",main="Environment discrimination capacity")
arrows(bp,z$score_ci_low,bp,z$score_ci_high,angle=90,code=3,length=.03,col="#333333"); abline(h=c(15,35,60),lty=3,col="#555555")
legend("topleft",names(class_col),fill=class_col,bty="n",ncol=2,cex=.8); dev.off()

png(file.path(fig,"图7_候选稳健性分类.png"),width=1400,height=900,res=170)
M<-t(as.matrix(rob_sum[c("robust","provisional","weak","not_screenable")])); colnames(M)<-rob_sum$nursery
barplot(M,col=c("#1b7837","#a6dba0","#fddbc7","#bdbdbd"),border=NA,ylab="Number of materials",main="Candidate robustness across six rules",legend.text=rownames(M),args.legend=list(x="topright",bty="n",cex=.8)); dev.off()

png(file.path(fig,"图8_区分力环境排名变化.png"),width=1600,height=1000,res=170)
par(mar=c(6,5,3,1)); R<-rbind(rank_pair_df$mean_abs_percentile_change,rank_pair_df$top_to_bottom_prop,1-rank_pair_df$spearman_rho)
colnames(R)<-rank_pair_df$nursery; rownames(R)<-c("Mean absolute percentile change","Top-bottom reversal proportion","1 - Spearman rho")
barplot(R,beside=TRUE,col=c("#fdae61","#d7191c","#2c7bb6"),border=NA,ylim=c(0,1),ylab="Instability metric (0-1)",main="Rank instability between discriminating sowing environments",legend.text=rownames(R),args.legend=list(x="topleft",bty="n",cex=.8))
dev.off()

png(file.path(fig,"图9_材料层面播期响应.png"),width=1900,height=1500,res=170)
groups<-split(sow_pair_df,interaction(sow_pair_df$nursery,sow_pair_df$location,sow_pair_df$year,drop=TRUE)); par(mfrow=c(2,4),mar=c(4,4,3,1))
resp_cols<-c(both_low="#1b7837",second_sowing_increase="#d73027",second_sowing_decrease="#4575b4",both_high="#7f0000",moderate_or_mixed="#bdbdbd")
for(d in groups){plot(d$first,d$second,pch=16,cex=.45,col=adjustcolor(resp_cols[d$response_class],.55),xlim=c(0,100),ylim=c(0,100),xlab="First sowing (%)",ylab="Second sowing (%)",main=paste(d$nursery[1],d$location[1],d$year[1]));abline(0,1,lty=2,col="#555555")};dev.off()

# Evidence-led phase 2 report.
envclass <- as.data.frame(table(env_disc$score_class)); names(envclass)<-c("class","n")
robtxt <- apply(rob_sum,1,function(x)sprintf("- %s：可评价%s个；稳健%s个，暂定%s个，弱证据%s个，不可评价%s个。",x["nursery"],x["eligible"],x["robust"],x["provisional"],x["weak"],x["not_screenable"]))
ranktxt <- apply(rank_pair_df,1,function(x)sprintf("- %s：%s vs %s，Spearman ρ=%s，平均绝对百分位变化=%s，前20%%到后20%%翻转%s个。",x["nursery"],x["environment1"],x["environment2"],fmt(as.numeric(x["spearman_rho"])),fmt(as.numeric(x["mean_abs_percentile_change"])),x["top_to_bottom_n"]))
cross_long <- aggregate(material_id~nursery+sowing_response_pattern,sow_mat_df,length)
sowtxt <- apply(cross_long,1,function(x)sprintf("- %s，%s：%s个材料。",x["nursery"],x["sowing_response_pattern"],x["material_id"]))
rob_sow <- aggregate(material_id~nursery+sowing_response_pattern,sow_mat_df[sow_mat_df$robustness_class=="robust",],length)
rob_sow_txt <- apply(rob_sow,1,function(x)sprintf("- %s稳健候选中，%s为%s类。",x["nursery"],x["material_id"],x["sowing_response_pattern"]))
report <- c("# 第二阶段：环境区分力、候选稳健性、排名翻转与播期响应","",sprintf("分析日期：%s  ",Sys.Date()),"数据集DOI：10.71682/10549333","",
"## 1. 环境区分力","",
sprintf("区分力评分由非零率、IQR、90–10百分位跨度和五级病害分布熵等权组成。19个圃×环境组合中：%s。",paste(sprintf("%s %s个",envclass$class,envclass$n),collapse="；")),
sprintf("三套判定规则中至少两套同意的共识区分力环境共有%s个。评分及1000次材料bootstrap区间见结果表15。该分数是探索性筛选工具，不是经过外部验证的量表。",sum(env_disc$consensus_discriminating)),"",
"## 2. 候选稳健性","","六套规则分别考察严格原始阈值、放宽阈值、环境内相对排名、全部环境最差值、圃内综合得分前20%以及无严重失败。入选≥4套定义为稳健，2–3套为暂定：",robtxt,"",
"41SAWSN缺少共识区分力环境，因此不把其低病害材料直接命名为稳定抗病候选。其他育种圃的稳健材料仍是后续验证优先级，不是遗传抗性证明。六套规则部分共享阈值信息，其一致性表示对规则扰动较稳健，不等于六份独立证据。","",
"## 3. 排名翻转","",ranktxt,"",
"排名变化只在共识区分力环境之间计算。百分位得分越高表示环境内病害越低。前/后20%集合按原始病害指数的20%和80%分位点划分，并保留分位点上的全部并列材料，因此集合可能超过精确的20%；大量并列值仍会降低排名分辨率。材料级环境特异类型见结果表19。","",
"## 4. 材料层面播期响应","",sowtxt,"","将候选稳健性与播期响应交叉后：",rob_sow_txt,"",
"播期结果以同一材料、地点和年份严格配对。`second_sowing_increase`要求第一播期≤5%而第二播期>10%；反向定义用于`second_sowing_decrease`。这些是表型响应类型，不作播期因果推断。","",
"## 对文章逻辑的贡献","","结果支持将文章主线设为：先评价筛选环境是否产生足够病害压力和材料区分度，再检验材料排名跨环境与跨播期的稳定性，最后识别对多套规则稳健的低病害候选。这样可以把“低病害”与“低压力下的逃病”明确区分。","",
"## 结果文件","","- `15_环境区分力综合评价.csv`：评分、bootstrap区间和三规则共识；","- `16_候选材料六规则稳健性矩阵.csv`：每个材料的六套选择结果；","- `18_环境间排名翻转汇总.csv`和`19_材料排名稳定性与环境特异性.csv`；","- `20_材料地点层面播期响应.csv`和`21_材料播期响应综合分类.csv`。")
writeLines(report,file.path(repd,"第二阶段四项深化分析报告.md"),useBytes=TRUE)
sink(file.path(logd,"第二阶段_R_sessionInfo.txt"));sessionInfo();sink()
cat("Phase 2 analyses completed\n")
