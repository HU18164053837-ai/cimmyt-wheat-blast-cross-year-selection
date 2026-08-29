options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260811)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_root <- file.path(root, "CIMMYT数据集")
tab_dir <- file.path(out_root, "结果表")
fig_dir <- file.path(out_root, "图表")
rep_dir <- file.path(out_root, "报告")
log_dir <- file.path(out_root, "运行记录")
for (d in c(tab_dir, fig_dir, rep_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

input_all <- file.path(root, "外部验证数据", "CIMMYT_2024", "处理后数据", "cimmyt_wheat_blast_long_all.csv")
input_obs <- file.path(root, "外部验证数据", "CIMMYT_2024", "处理后数据", "cimmyt_wheat_blast_long_observed.csv")
stopifnot(file.exists(input_all), file.exists(input_obs))

all_dat <- read.csv(input_all, check.names = FALSE, na.strings = c(""), fileEncoding = "UTF-8-BOM")
obs <- read.csv(input_obs, check.names = FALSE, na.strings = c(""), fileEncoding = "UTF-8-BOM")
names(all_dat)[names(all_dat) == "Ent"] <- "Entry"
names(obs)[names(obs) == "Ent"] <- "Entry"
all_dat$disease_index_pct <- suppressWarnings(as.numeric(all_dat$disease_index_pct))
obs$disease_index_pct <- suppressWarnings(as.numeric(obs$disease_index_pct))
all_dat$Entry <- as.character(all_dat$Entry)
obs$Entry <- as.character(obs$Entry)
all_dat$GID <- as.character(all_dat$GID)
obs$GID <- as.character(obs$GID)
all_dat$material_id <- paste(all_dat$nursery, all_dat$Entry, sep = "::")
obs$material_id <- paste(obs$nursery, obs$Entry, sep = "::")
all_dat$env_id <- paste(all_dat$location, all_dat$year, all_dat$sowing, sep = "_")
obs$env_id <- paste(obs$location, obs$year, obs$sowing, sep = "_")

stopifnot(all(obs$disease_index_pct >= 0 & obs$disease_index_pct <= 100))
nurseries <- sort(unique(obs$nursery))

qfun <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))
safe_num <- function(x) ifelse(is.finite(x), x, NA_real_)
write_csv <- function(x, name) write.csv(x, file.path(tab_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")

# 1. Data overview and environment pressure
overview <- data.frame(
  nursery = nurseries,
  material_rows = sapply(nurseries, function(z) length(unique(all_dat$material_id[all_dat$nursery == z]))),
  environments = sapply(nurseries, function(z) length(unique(all_dat$env_id[all_dat$nursery == z]))),
  potential_cells = sapply(nurseries, function(z) sum(all_dat$nursery == z)),
  observed_cells = sapply(nurseries, function(z) sum(obs$nursery == z)),
  missing_cells = sapply(nurseries, function(z) sum(all_dat$nursery == z & is.na(all_dat$disease_index_pct)))
)
overview$missing_rate <- overview$missing_cells / overview$potential_cells
write_csv(overview, "01_数据概况.csv")

env_keys <- unique(all_dat[c("nursery", "environment", "env_id", "location", "year", "sowing")])
env_list <- lapply(seq_len(nrow(env_keys)), function(i) {
  k <- env_keys[i, ]
  a <- all_dat[all_dat$nursery == k$nursery & all_dat$environment == k$environment, ]
  x <- a$disease_index_pct
  y <- x[!is.na(x)]
  data.frame(k, n_materials = nrow(a), n_observed = length(y), n_missing = sum(is.na(x)),
             missing_rate = mean(is.na(x)), zero_prop = mean(y == 0), ge20_prop = mean(y >= 20),
             mean = mean(y), sd = sd(y), median = median(y), q1 = qfun(y, .25), q3 = qfun(y, .75), max = max(y))
})
env_summary <- do.call(rbind, env_list)
env_summary <- env_summary[order(env_summary$nursery, -env_summary$mean), ]
env_summary$informative_environment <- env_summary$zero_prop < .80 & env_summary$q3 > 0
row.names(env_summary) <- NULL
write_csv(env_summary, "02_环境病害压力汇总.csv")
informative_keys <- paste(env_summary$nursery[env_summary$informative_environment],
                          env_summary$environment[env_summary$informative_environment], sep="::")

# Within-environment relative resistance score, ties retained.
obs$relative_resistance <- NA_real_
for (idx in split(seq_len(nrow(obs)), interaction(obs$nursery, obs$environment, drop = TRUE))) {
  x <- obs$disease_index_pct[idx]
  if (length(x) <= 1L) obs$relative_resistance[idx] <- 1 else {
    obs$relative_resistance[idx] <- 1 - (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  }
}

# 2. Material summaries and transparent candidate rule
mat_meta <- unique(all_dat[c("nursery", "Entry", "GID", "Cross", "Sel_Hist", "material_id")])
mat_rows <- lapply(split(obs, obs$material_id), function(d) {
  x <- d$disease_index_pct
  inf <- d[paste(d$nursery,d$environment,sep="::") %in% informative_keys,]
  xi <- inf$disease_index_pct
  n_total <- length(unique(all_dat$environment[all_dat$material_id == d$material_id[1]]))
  data.frame(nursery = d$nursery[1], Entry = d$Entry[1], GID = d$GID[1], Cross = d$Cross[1],
             Sel_Hist = d$Sel_Hist[1], material_id = d$material_id[1], n_observed = length(x),
             n_total_env = n_total, complete_rate = length(x) / n_total, mean_disease = mean(x),
             median_disease = median(x), max_disease = max(x), sd_disease = if(length(x)>1) sd(x) else NA,
             range_disease = max(x)-min(x), low5_prop = mean(x <= 5), ge20_prop = mean(x >= 20),
             mean_relative_resistance = mean(d$relative_resistance),
             worst_relative_resistance = min(d$relative_resistance), n_informative_env=length(xi),
             informative_mean=if(length(xi)) mean(xi) else NA,
             informative_max=if(length(xi)) max(xi) else NA,
             informative_low5_prop=if(length(xi)) mean(xi<=5) else NA)
})
material_summary <- do.call(rbind, mat_rows)
material_summary$stable_low_candidate <- with(material_summary,
  complete_rate >= .80 & n_informative_env >= 2 & informative_mean <= 5 &
    informative_max <= 10 & informative_low5_prop >= .80)
material_summary$composite_score <- 100 * (0.50 * material_summary$mean_relative_resistance +
  0.25 * material_summary$worst_relative_resistance + 0.25 * material_summary$low5_prop)
material_summary <- material_summary[order(material_summary$nursery, -material_summary$stable_low_candidate,
                                             -material_summary$composite_score, material_summary$mean_disease), ]
row.names(material_summary) <- NULL
write_csv(material_summary, "03_材料跨环境表现与候选筛选.csv")

candidate_counts <- do.call(rbind, lapply(split(material_summary, material_summary$nursery), function(d) {
  ni <- sum(env_summary$nursery==d$nursery[1] & env_summary$informative_environment)
  data.frame(nursery=d$nursery[1], n_materials=nrow(d), complete80=sum(d$complete_rate>=.8),
             informative_environments=ni, screenable=ni>=2,
             stable_low_candidates=sum(d$stable_low_candidate), candidate_rate=mean(d$stable_low_candidate))
}))
write_csv(candidate_counts, "04_候选材料数量汇总.csv")

top_candidates <- do.call(rbind, lapply(split(material_summary, material_summary$nursery), function(d) {
  d <- d[order(-d$stable_low_candidate, -d$composite_score, d$mean_disease), ]
  head(d, 20)
}))
write_csv(top_candidates, "05_各育种圃前20名候选材料.csv")

# 3. Environment correlations
cor_rows <- list(); cix <- 0L
for (nur in nurseries) {
  d <- obs[obs$nursery == nur, c("material_id", "env_id", "disease_index_pct")]
  envs <- sort(unique(d$env_id))
  for (i in seq_along(envs)) for (j in i:length(envs)) {
    a <- d[d$env_id == envs[i], c("material_id", "disease_index_pct")]
    b <- d[d$env_id == envs[j], c("material_id", "disease_index_pct")]
    names(a)[2] <- "x"; names(b)[2] <- "y"; m <- merge(a,b,by="material_id")
    rho <- if(nrow(m)>=3 && sd(m$x)>0 && sd(m$y)>0) cor(m$x,m$y,method="spearman") else NA
    cix <- cix+1L; cor_rows[[cix]] <- data.frame(nursery=nur, environment1=envs[i], environment2=envs[j], n_pairs=nrow(m), spearman_rho=rho)
  }
}
env_cor <- do.call(rbind, cor_rows)
write_csv(env_cor, "06_环境间Spearman相关.csv")

# 4. Paired sowing effects with deterministic paired bootstrap CIs
boot_mean_ci <- function(x, B=2000L) {
  if(length(x)<2) return(c(NA,NA))
  bs <- replicate(B, mean(sample(x, length(x), replace=TRUE)))
  as.numeric(quantile(bs, c(.025,.975), names=FALSE))
}
pairs_key <- unique(obs[c("nursery","location","year")])
sow_rows <- list(); six <- 0L
paired_values <- list()
for (i in seq_len(nrow(pairs_key))) {
  k <- pairs_key[i,]
  d <- obs[obs$nursery==k$nursery & obs$location==k$location & obs$year==k$year,]
  a <- d[d$sowing=="1st", c("material_id","disease_index_pct")]
  b <- d[d$sowing=="2nd", c("material_id","disease_index_pct")]
  if(!nrow(a) || !nrow(b)) next
  names(a)[2] <- "first"; names(b)[2] <- "second"; m <- merge(a,b,by="material_id")
  dif <- m$second-m$first; ci <- boot_mean_ci(dif)
  wt <- tryCatch(wilcox.test(m$second,m$first,paired=TRUE,exact=FALSE), error=function(e) NULL)
  six <- six+1L
  sow_rows[[six]] <- data.frame(nursery=k$nursery, location=k$location, year=k$year, n_pairs=nrow(m),
    first_mean=mean(m$first), second_mean=mean(m$second), mean_difference=mean(dif),
    mean_diff_ci_low=ci[1], mean_diff_ci_high=ci[2], median_difference=median(dif),
    prop_second_higher=mean(dif>0), prop_equal=mean(dif==0), wilcoxon_p=if(is.null(wt)) NA else wt$p.value)
  paired_values[[six]] <- data.frame(nursery=k$nursery,location=k$location,year=k$year,material_id=m$material_id,difference=dif)
}
sowing_summary <- do.call(rbind,sow_rows)
sowing_summary$wilcoxon_p_bh <- p.adjust(sowing_summary$wilcoxon_p, method="BH")
write_csv(sowing_summary, "07_播期材料内配对比较.csv")
paired_long <- do.call(rbind,paired_values)
write_csv(paired_long, "08_播期配对差明细.csv")

# 5. Finlay-Wilkinson descriptive stability
fw_rows <- list(); fix <- 0L
for(nur in nurseries) {
  d <- obs[obs$nursery==nur,]
  env_means <- aggregate(disease_index_pct~env_id,d,mean)
  names(env_means)[2] <- "env_mean"
  d <- merge(d,env_means,by="env_id")
  for(mid in unique(d$material_id)) {
    z <- d[d$material_id==mid,]
    if(nrow(z)>=3 && sd(z$env_mean)>0) {
      fit <- lm(disease_index_pct~env_mean,z)
      fix <- fix+1L
      fw_rows[[fix]] <- data.frame(nursery=nur,material_id=mid,Entry=z$Entry[1],GID=z$GID[1],n_env=nrow(z),
        intercept=coef(fit)[1],environment_slope=coef(fit)[2],residual_mse=mean(residuals(fit)^2),r_squared=summary(fit)$r.squared)
    }
  }
}
fw <- do.call(rbind,fw_rows)
write_csv(fw,"09_Finlay_Wilkinson稳定性.csv")

# 6. Exploratory AMMI SVD on complete material matrices
ammi_gen <- list(); ammi_env <- list(); ammi_var <- list(); aix <- 0L
for(nur in nurseries) {
  d <- obs[obs$nursery==nur,c("material_id","env_id","disease_index_pct")]
  mids <- sort(unique(d$material_id)); envs <- sort(unique(d$env_id))
  M <- matrix(NA_real_, length(mids), length(envs), dimnames=list(mids,envs))
  M[cbind(match(d$material_id,mids),match(d$env_id,envs))] <- d$disease_index_pct
  M <- M[complete.cases(M),,drop=FALSE]
  if(nrow(M)<3 || ncol(M)<2) next
  grand <- mean(M); R <- sweep(M,1,rowMeans(M),"-"); R <- sweep(R,2,colMeans(M),"-"); R <- R + grand
  sv <- svd(R); vv <- sv$d^2/sum(sv$d^2)
  gsc <- sv$u %*% diag(sqrt(sv$d), nrow=length(sv$d)); esc <- sv$v %*% diag(sqrt(sv$d), nrow=length(sv$d))
  aix <- aix+1L
  ammi_gen[[aix]] <- data.frame(nursery=nur,material_id=rownames(M),PC1=gsc[,1],PC2=if(ncol(gsc)>=2) gsc[,2] else NA)
  ammi_env[[aix]] <- data.frame(nursery=nur,environment=colnames(M),PC1=esc[,1],PC2=if(ncol(esc)>=2) esc[,2] else NA)
  ammi_var[[aix]] <- data.frame(nursery=nur,n_complete_materials=nrow(M),n_environments=ncol(M),PC=seq_along(vv),explained_interaction=vv)
}
ammi_g <- do.call(rbind,ammi_gen); ammi_e <- do.call(rbind,ammi_env); ammi_v <- do.call(rbind,ammi_var)
write_csv(ammi_g,"10_AMMI材料得分.csv"); write_csv(ammi_e,"11_AMMI环境得分.csv"); write_csv(ammi_v,"12_AMMI解释比例.csv")

# 7. Shared GIDs and repeated local checks
gid_nur <- unique(mat_meta[c("GID","nursery")])
shared_ids <- names(which(table(gid_nur$GID)>=2))
shared_rows <- lapply(shared_ids,function(g) {
  d <- obs[obs$GID==g,]
  byn <- aggregate(disease_index_pct~nursery,d,mean)
  data.frame(GID=g,n_nurseries=length(unique(d$nursery)),nurseries=paste(sort(unique(d$nursery)),collapse=";"),
             n_observations=nrow(d),overall_mean=mean(d$disease_index_pct),nursery_mean_min=min(byn$disease_index_pct),
             nursery_mean_max=max(byn$disease_index_pct),nursery_mean_range=max(byn$disease_index_pct)-min(byn$disease_index_pct))
})
shared_summary <- do.call(rbind,shared_rows)
write_csv(shared_summary,"13_跨圃共享GID表现.csv")

dup_keys <- names(which(table(paste(mat_meta$nursery,mat_meta$GID,sep="::"))>1))
dup_detail <- list(); dix <- 0L
for(key in dup_keys) {
  sp <- strsplit(key,"::",fixed=TRUE)[[1]]; d <- obs[obs$nursery==sp[1]&obs$GID==sp[2],]
  for(e in unique(d$env_id)) {
    z <- d[d$env_id==e,]
    if(nrow(z)>=2) {
      vals <- z$disease_index_pct
      dix <- dix+1L; dup_detail[[dix]] <- data.frame(nursery=sp[1],GID=sp[2],environment=e,n_entries=nrow(z),
        min=min(vals),max=max(vals),absolute_range=max(vals)-min(vals),all_local_check=all(grepl("LOCAL CHECK",z$Cross,fixed=TRUE)|grepl("LOCAL CHECK",z$Sel_Hist,fixed=TRUE)))
    }
  }
}
dup_check <- if(length(dup_detail)) do.call(rbind,dup_detail) else data.frame()
write_csv(dup_check,"14_重复GID对照一致性.csv")

# Figures: base R only, English labels for portable rendering
png(file.path(fig_dir,"图1_环境病害压力.png"),width=1800,height=1200,res=170)
par(mar=c(8,5,3,1)); ord <- order(env_summary$mean)
loc3 <- substr(env_summary$location,1,3); labs <- paste0(env_summary$nursery,"\n",loc3,substr(env_summary$year,3,4),"-S",substr(env_summary$sowing,1,1))
bpraw <- lapply(seq_len(nrow(env_summary))[ord],function(i) obs$disease_index_pct[obs$nursery==env_summary$nursery[i]&obs$env_id==env_summary$env_id[i]])
bpdat <- lapply(bpraw,sqrt)
boxplot(bpdat,names=labs[ord],las=2,outline=FALSE,col="#9ecae1",border="#2b6f8a",yaxt="n",ylab="Wheat blast index (%) [square-root display]",main="Disease pressure across nursery-environments",ylim=c(0,10))
axis(2,at=sqrt(c(0,1,5,10,20,50,100)),labels=c(0,1,5,10,20,50,100))
stripchart(bpdat,vertical=TRUE,method="jitter",pch=16,cex=.25,col=adjustcolor("#08306b",.35),add=TRUE)
dev.off()

png(file.path(fig_dir,"图2_材料平均值与最差环境.png"),width=1800,height=1400,res=170)
par(mfrow=c(2,2),mar=c(5,5,3,1))
for(nur in nurseries) {
 d<-material_summary[material_summary$nursery==nur,]; cols<-ifelse(d$stable_low_candidate,"#17823b",adjustcolor("#666666",.45))
 plot(d$mean_disease,d$max_disease,pch=16,col=cols,xlab="Mean disease index (%)",ylab="Maximum disease index (%)",main=nur)
 abline(v=5,h=10,lty=2,col="#b2182b"); legend("topleft",c("Rule-based candidate","Other"),pch=16,col=c("#17823b","#999999"),bty="n",cex=.8)
}
dev.off()

png(file.path(fig_dir,"图3_播期配对差.png"),width=1500,height=950,res=170)
par(mar=c(8,5,3,1)); paired_long$group <- paste(paired_long$nursery,paired_long$location,paired_long$year,sep=" | ")
boxplot(difference~group,paired_long,las=2,col="#fdae6b",border="#a63603",ylab="Second - first sowing disease index (%)",main="Within-material sowing-date differences",outline=FALSE)
abline(h=0,lty=2,col="#444444"); stripchart(difference~group,paired_long,vertical=TRUE,method="jitter",add=TRUE,pch=16,cex=.22,col=adjustcolor("#7f2704",.3))
dev.off()

png(file.path(fig_dir,"图4_AMMI环境互作双标图.png"),width=1800,height=1400,res=170)
par(mfrow=c(2,2),mar=c(5,5,3,1))
for(nur in nurseries) {
 g<-ammi_g[ammi_g$nursery==nur,]; e<-ammi_e[ammi_e$nursery==nur,]; vr<-ammi_v[ammi_v$nursery==nur,]
 xr<-range(c(g$PC1,e$PC1),na.rm=TRUE); yr<-range(c(g$PC2,e$PC2),na.rm=TRUE); xpad=max(diff(xr)*.12,.5); ypad=max(diff(yr)*.12,.5); xlim=xr+c(-xpad,xpad); ylim=yr+c(-ypad,ypad)
 plot(g$PC1,g$PC2,pch=16,cex=.45,col=adjustcolor("#636363",.45),xlim=xlim,ylim=ylim,xlab=sprintf("PC1 (%.1f%%)",100*vr$explained_interaction[1]),ylab=sprintf("PC2 (%.1f%%)",100*vr$explained_interaction[2]),main=nur)
 shorte <- gsub("Quirusillas","Qui",gsub("Jashore","Jas",gsub("Okinawa","Oki",e$environment)))
 abline(h=0,v=0,col="#bdbdbd"); points(e$PC1,e$PC2,pch=17,col="#cb181d",cex=1.1); text(e$PC1,e$PC2,labels=shorte,pos=3,cex=.52,col="#99000d",xpd=NA)
}
dev.off()

png(file.path(fig_dir,"图5_环境相关性热图.png"),width=1800,height=1400,res=170)
par(mfrow=c(2,2),mar=c(8,8,3,2))
for(nur in nurseries) {
 d<-env_cor[env_cor$nursery==nur,]; es<-sort(unique(c(d$environment1,d$environment2))); M<-matrix(NA,length(es),length(es),dimnames=list(es,es))
 for(i in seq_len(nrow(d))) { M[d$environment1[i],d$environment2[i]]<-d$spearman_rho[i]; M[d$environment2[i],d$environment1[i]]<-d$spearman_rho[i] }
 image(seq_along(es),seq_along(es),t(M[nrow(M):1,,drop=FALSE]),zlim=c(-1,1),col=colorRampPalette(c("#2166ac","white","#b2182b"))(101),axes=FALSE,xlab="",ylab="",main=nur)
 axis(1,seq_along(es),es,las=2,cex.axis=.55); axis(2,seq_along(es),rev(es),las=2,cex.axis=.55)
}
dev.off()

# Compact evidence-led report
fmt <- function(x,d=2) ifelse(is.na(x),"NA",format(round(x,d),nsmall=d,trim=TRUE))
highest <- env_summary[which.max(env_summary$mean),]
lowest <- env_summary[which.min(env_summary$mean),]
pc1 <- ammi_v[ammi_v$PC==1,c("nursery","explained_interaction")]
sow_text <- apply(sowing_summary,1,function(x) sprintf("- %s，%s（%s）：配对数%s，第二播期−第一播期均值差%s个百分点（bootstrap 95%%区间%s至%s），中位数差%s，BH校正P=%s。",x["nursery"],x["location"],x["year"],x["n_pairs"],fmt(as.numeric(x["mean_difference"])),fmt(as.numeric(x["mean_diff_ci_low"])),fmt(as.numeric(x["mean_diff_ci_high"])),fmt(as.numeric(x["median_difference"])),format.pval(as.numeric(x["wilcoxon_p_bh"]),digits=3,eps=1e-4)))
cand_text <- apply(candidate_counts,1,function(x) if(x["screenable"]=="TRUE") sprintf("- %s：有%s个有区分力环境；%s/%s个材料满足预设候选规则（%s%%）。",x["nursery"],x["informative_environments"],x["stable_low_candidates"],x["n_materials"],fmt(100*as.numeric(x["candidate_rate"]),1)) else sprintf("- %s：仅%s个有区分力环境，当前不可作稳定候选筛选；保留材料排序供后续试验参考。",x["nursery"],x["informative_environments"]))
report <- c(
"# CIMMYT小麦穗瘟病多环境表型第一阶段分析报告",
"",sprintf("分析日期：%s  ",Sys.Date()),"数据集DOI：10.71682/10549333","",
"## 核心结论","",
sprintf("本阶段分析覆盖4个育种圃、%s个材料行、%s个圃×环境组合和%s个非缺失病害指数。分析是探索性多环境表型评价，不是气象模型外部验证。",sum(overview$material_rows),sum(overview$environments),sum(overview$observed_cells)),
sprintf("最高平均病害压力出现在%s的%s（均值%s%%，中位数%s%%）；最低出现在%s的%s（均值%s%%）。环境间压力高度不均衡，多个环境呈大量零值，因此单一低压环境中的零病害不能单独证明稳定抗性。",highest$nursery,highest$env_id,fmt(highest$mean),fmt(highest$median),lowest$nursery,lowest$env_id,fmt(lowest$mean)),
"", "有区分力环境定义为零值比例<80%且上四分位数>0。候选需覆盖至少2个有区分力环境，并满足完整率≥80%、有区分力环境均值≤5%、最大值≤10且至少80%观测≤5：",cand_text,"",
"## 播期响应","",sow_text,"",
"播期比较严格限定为同一育种圃、材料、地点和年份的配对结果。正差表示第二播期病害更高。由于播期并非随机处理，结果描述关联而非因果效应。","",
"## 材料×环境互作","",
paste0("AMMI完整矩阵中PC1对交互变异的解释比例分别为：",paste(sprintf("%s %.1f%%",pc1$nursery,100*pc1$explained_interaction),collapse="；"),"。AMMI仅用于展示交互结构；公开数据缺少区组与小区重复，不能由此估计田间误差或作正式显著性检验。"),
"环境间Spearman相关、Finlay–Wilkinson斜率和AMMI得分均已完整输出。环境中大量并列零值会压缩相关和排名信息，解释时应同时参考环境病害压力。","",
"## 数据质量与身份核对","",
sprintf("共识别%s个跨至少两个育种圃出现的GID；它们用于身份与表型一致性描述，不直接合并为同一试验单位。重复GID对照在共享环境中的一致性见结果表14。",length(shared_ids)),
"15HZAN的环境列标记为2023，而数据集标题称2024周期。本分析忠实保留列名年份；在CIMMYT确认前，不将其改写为2024。","",
"## 候选材料的使用方式","",
"结果表03包含所有材料指标和规则标记，结果表05给出各育种圃前20名。候选规则是透明的分析阈值，不是官方抗性分级。优先关注同时具备低均值、低最大值、高完整率和较高环境内相对得分的材料，并在独立重复试验中验证。","",
"## 主要限制","",
"- 缺少明确区组、重复、小区、接种方式和调查时间信息；","- 环境数有限且病害指数零膨胀；","- 4个育种圃材料集合不同，主体结论必须在圃内解释；","- 分析为事后探索，P值仅作辅助并进行了BH校正；","- 未经CIMMYT书面许可，不应重新分发原始数据或可完整还原原始数据的长表。","",
"## 文件导航","",
"- `结果表/02_环境病害压力汇总.csv`：环境压力；","- `结果表/03_材料跨环境表现与候选筛选.csv`：全体材料与候选规则；","- `结果表/07_播期材料内配对比较.csv`：播期效应；","- `结果表/09_Finlay_Wilkinson稳定性.csv`：稳定性斜率；","- `结果表/10_AMMI材料得分.csv`至`12_AMMI解释比例.csv`：交互结构；","- `结果表/13_跨圃共享GID表现.csv`与`14_重复GID对照一致性.csv`：身份和重复对照检查。"
)
writeLines(report,file.path(rep_dir,"第一阶段分析报告.md"),useBytes=TRUE)

# Reproducibility records
hashes <- tools::md5sum(c(input_all,input_obs))
writeLines(c(sprintf("analysis_time=%s",Sys.time()),sprintf("input=%s md5=%s",basename(names(hashes)),unname(hashes))),file.path(log_dir,"输入文件校验.txt"),useBytes=TRUE)
sink(file.path(log_dir,"R_sessionInfo.txt")); sessionInfo(); sink()
cat("Phase 1 analysis completed. Outputs:", out_root, "\n")
