############################### UVMR-ALT/RETINA #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(BWMR)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_RETINA")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALT_RETINA/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_RETINOPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_RETINA/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_RETINA/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_RETINA/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_RETINA/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_RETINA/MR_result.csv")
pdf(file = "./result/UVMR_ALT_RETINA/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_RETINA/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_RETINA/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_RETINA/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_RETINA/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_RETINA/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])


############################### UVMR-ALT/RETINA1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_RETINA1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALT_RETINA1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_RETINA1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_RETINA1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_RETINA1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_RETINA1/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA1/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_RETINA1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_RETINA1/MR_result.csv")
pdf(file = "./result/UVMR_ALT_RETINA1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_RETINA1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_RETINA1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_RETINA1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_RETINA1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_RETINA1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])

############################### UVMR-ALT/RETINA2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_RETINA2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALT_RETINA2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_RETINA2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_RETINA2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA2/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs498475", "rs1126506", "rs4084164", "rs2642438", "rs4135250")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_RETINA2/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_RETINA2/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_RETINA2/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_RETINA2/MR_result.csv")
pdf(file = "./result/UVMR_ALT_RETINA2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_RETINA2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_RETINA2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_RETINA2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_RETINA2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_RETINA2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])

############################### UVMR-ALT/NEPHRO #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_NEPHRO")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000)
write.csv(exposure, file = "./result/UVMR_ALT_NEPHRO/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_NEPHROPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_NEPHRO/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_NEPHRO/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_NEPHRO/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_NEPHRO/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs34082903", "rs4841133", "rs686250", "rs429358")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_NEPHRO/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_NEPHRO/MR_result.csv")
pdf(file = "./result/UVMR_ALT_NEPHRO/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_NEPHRO/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_NEPHRO/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_NEPHRO/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_NEPHRO/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_NEPHRO/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])

############################### UVMR-ALT/NEPHRO1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_NEPHRO1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure$SNP, file = "./result/UVMR_ALT_NEPHRO1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_NEPHRO1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_NEPHRO1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_NEPHRO1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_NEPHRO1/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO1/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_NEPHRO1/MR_result.csv")
pdf(file = "./result/UVMR_ALT_NEPHRO1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_NEPHRO1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_NEPHRO1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_NEPHRO1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_NEPHRO1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_NEPHRO1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])

############################### UVMR-ALT/NEPHRO2 #####################################
#打开R包``````````````````````
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALT_NEPHRO2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30620_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALT_NEPHRO2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALT_NEPHRO2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALT_NEPHRO2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALT_NEPHRO2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALT_NEPHRO2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO2/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs35149321", "rs429358")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALT_NEPHRO2/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALT_NEPHRO2/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_ALT_NEPHRO2/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALT_NEPHRO2/MR_result.csv")
pdf(file = "./result/UVMR_ALT_NEPHRO2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALT_NEPHRO2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALT_NEPHRO2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALT_NEPHRO2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALT_NEPHRO2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALT_NEPHRO2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])

############################### UVMR-AST/RETINA #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_RETINA")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_RETINA/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_RETINOPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_RETINA/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_RETINA/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_RETINA/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs628401")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_RETINA/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_RETINA/MR_result.csv")
pdf(file = "./result/UVMR_AST_RETINA/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_RETINA/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_RETINA/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_RETINA/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_RETINA/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_RETINA/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-AST/RETINA1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_RETINA1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_RETINA1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_RETINA1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_RETINA1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_RETINA1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_RETINA1/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA1/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs738409", "rs12359178")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_RETINA1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_RETINA1/MR_result.csv")
pdf(file = "./result/UVMR_AST_RETINA1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_RETINA1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_RETINA1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_RETINA1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_RETINA1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_RETINA1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-AST/RETINA2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_RETINA2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_RETINA2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_RETINA2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_RETINA2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA2/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs201376697", "rs899729", "rs11714389", "rs57562692")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA2/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA2/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_RETINA2/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_RETINA2/MR_result.csv")
pdf(file = "./result/UVMR_AST_RETINA2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_RETINA2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_RETINA2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_RETINA2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_RETINA2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_RETINA2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-AST/NEPHRO #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_NEPHRO")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_NEPHRO/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_NEPHROPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_NEPHRO/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_NEPHRO/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_NEPHRO/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_NEPHRO/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs2367585")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_RETINA2/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_RETINA2/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_RETINA2/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_NEPHRO/MR_result.csv")
pdf(file = "./result/UVMR_AST_NEPHRO/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_NEPHRO/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_NEPHRO/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_NEPHRO/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_NEPHRO/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_NEPHRO/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-AST/NEPHRO1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_NEPHRO1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_NEPHRO1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_NEPHRO1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_NEPHRO1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_NEPHRO1/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs12359178")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_NEPHRO1/MR_result.csv")
pdf(file = "./result/UVMR_AST_NEPHRO1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_NEPHRO1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_NEPHRO1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_NEPHRO1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_NEPHRO1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_NEPHRO1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-AST/NEPHRO2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_AST_NEPHRO2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30650_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_AST_NEPHRO2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_AST_NEPHRO2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_AST_NEPHRO2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_AST_NEPHRO2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_AST_NEPHRO2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO2/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_AST_NEPHRO2/MR_result.csv")
pdf(file = "./result/UVMR_AST_NEPHRO2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_AST_NEPHRO2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_AST_NEPHRO2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_AST_NEPHRO2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_AST_NEPHRO2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_AST_NEPHRO2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=10.65, p=MR_result[1,9])

############################### UVMR-ALP/RETINA #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_RETINA")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_RETINA/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_RETINOPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_RETINA/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_RETINA/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_RETINA/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_RETINA/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_RETINA/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs2393791")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_RETINA/MR_result.csv")
pdf(file = "./result/UVMR_ALP_RETINA/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_RETINA/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_RETINA/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_RETINA/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_RETINA/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_RETINA/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-ALP/RETINA1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_RETINA1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_RETINA1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_RETINA1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_RETINA1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_RETINA1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_RETINA1/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_RETINA1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_RETINA1/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs6469788", "rs7855852", "rs114639894")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_RETINA1/MR_result.csv")
pdf(file = "./result/UVMR_ALP_RETINA1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_RETINA1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_RETINA1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_RETINA1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_RETINA1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_RETINA1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-ALP/RETINA2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_RETINA2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_RETINA2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_RETINA2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_RETINA2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_RETINA2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_RETINA2/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs603424")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_RETINA2/MR_result.csv")
pdf(file = "./result/UVMR_ALP_RETINA2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_RETINA2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_RETINA2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_RETINA2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_RETINA2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_RETINA2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-ALP/NEPHRO #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_NEPHRO")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_NEPHRO/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_NEPHROPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_NEPHRO/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_NEPHRO/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_NEPHRO/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_NEPHRO/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_NEPHRO/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_NEPHRO/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs56222534", "rs41282145", "rs2393791")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_NEPHRO/MR_result.csv")
pdf(file = "./result/UVMR_ALP_NEPHRO/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_NEPHRO/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_NEPHRO/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_NEPHRO/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_NEPHRO/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_NEPHRO/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-ALP/NEPHRO1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_NEPHRO1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_NEPHRO1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_NEPHRO1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_NEPHRO1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_NEPHRO1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_NEPHRO/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_NEPHRO1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_NEPHRO1/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_NEPHRO1/MR_result.csv")
pdf(file = "./result/UVMR_ALP_NEPHRO1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_NEPHRO1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_NEPHRO1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_NEPHRO1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_NEPHRO1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_NEPHRO1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-ALP/NEPHRO2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_ALP_NEPHRO2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30610_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_ALP_NEPHRO2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_ALP_NEPHRO2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_ALP_NEPHRO2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_ALP_NEPHRO2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_ALP_NEPHRO2/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_ALP_NEPHRO2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_ALP_NEPHRO2/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_ALP_NEPHRO2/MR_result.csv")
pdf(file = "./result/UVMR_ALP_NEPHRO2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_ALP_NEPHRO2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_ALP_NEPHRO2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_ALP_NEPHRO2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_ALP_NEPHRO2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_ALP_NEPHRO2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=26.38, p=MR_result[1,9])

############################### UVMR-GGT/RETINA #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_RETINA")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_RETINA/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_RETINOPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_RETINA/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_RETINA/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_RETINA/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_RETINA/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_RETINA/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs9738226")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_RETINA/MR_result.csv")
pdf(file = "./result/UVMR_GGT_RETINA/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_RETINA/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_RETINA/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_RETINA/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_RETINA/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_RETINA/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### UVMR-GGT/RETINA1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_RETINA1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_RETINA1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_RETINA1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_RETINA1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_RETINA1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_RETINA1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_RETINA1/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_RETINA1/MR_result.csv")
pdf(file = "./result/UVMR_GGT_RETINA1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_RETINA1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_RETINA1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_RETINA1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_RETINA1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_RETINA1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### UVMR-GGT/RETINA2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_RETINA2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_RETINA2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_RETINA2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_RETINA2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_RETINA2/MRPRESSO_outlier_result.csv")

mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_RETINA2/MR_result.csv")
pdf(file = "./result/UVMR_GGT_RETINA2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_RETINA2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_RETINA2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_RETINA2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_RETINA2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_RETINA2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### UVMR-GGT/NEPHRO #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_NEPHRO")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_NEPHRO/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_DM_NEPHROPATHY_EXMORE")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_NEPHRO/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_NEPHRO/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_NEPHRO/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_NEPHRO/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_NEPHRO/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs9738226")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_AST_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_AST_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_NEPHRO/MR_result.csv")
pdf(file = "./result/UVMR_GGT_NEPHRO/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_NEPHRO/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_NEPHRO/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_NEPHRO/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_NEPHRO/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_NEPHRO/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### UVMR-GGT/NEPHRO1 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_NEPHRO1")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_NEPHRO1/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM1REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_NEPHRO1/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_NEPHRO1/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_NEPHRO1/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs9738226", "rs3859862")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_GGT_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_NEPHRO1/MR_result.csv")
pdf(file = "./result/UVMR_GGT_NEPHRO1/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_NEPHRO1/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_NEPHRO1/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_NEPHRO1/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_NEPHRO1/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_NEPHRO1/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### UVMR-GGT/NEPHRO2 #####################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/UVMR_GGT_NEPHRO2")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ukb-d-30730_raw", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/UVMR_GGT_NEPHRO2/exposure_raw.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2REN")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/UVMR_GGT_NEPHRO2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/UVMR_GGT_NEPHRO2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/UVMR_GGT_NEPHRO2/data_harmonised.csv")

#去除与混淆因素相关的SNP
snp_excluded = read.table("./result/UVMR_GGT_RETINA/SNP_excluded.txt")
snp_excluded = snp_excluded$V1
data = data[!(data$SNP %in% snp_excluded),]
snp_excluded

#孟德尔随机化
#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_NEPHRO2/MRPRESSO_outlier.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_NEPHRO2/MRPRESSO_outlier_result.csv")

data = data[!(data$SNP %in% c("rs9738226", "rs3859862")),]
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso
save(mr_presso, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_final.Rdata")
write.csv(mr_presso[[1]]$`Main MR results`, file = "./result/UVMR_GGT_NEPHRO1/MRPRESSO_final_result.csv")
write.csv(data, file = "./result/UVMR_GGT_NEPHRO1/data_harmonised_final.csv")

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/UVMR_GGT_NEPHRO2/MR_result.csv")
pdf(file = "./result/UVMR_GGT_NEPHRO2/scatterPlot.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/UVMR_GGT_NEPHRO2/MR_heterogeneity.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/UVMR_GGT_NEPHRO2/MR_pleiotropy.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/UVMR_GGT_NEPHRO2/funnelPlot.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/UVMR_GGT_NEPHRO2/leaveoneoutPlot.pdf", width = 8, height = 16)
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res
save(BWMR_res, file = "./result/UVMR_GGT_NEPHRO2/BWMR_res.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=41.87, p=MR_result[1,9])

############################### MVMR-ALT-HbA1c-SBP/RETINA2 #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/MVMR_ALT-HbA1c-SBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30620_raw", "ebi-a-GCST90002244", "ieu-b-38")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_ALT-HbA1c-SBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_ALT-HbA1c-SBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs961764", "rs1012089", "rs12321", "rs13395911", "rs1821002",
                                        "rs1870735", "rs2024385", "rs2954027", "rs3802517", "rs3828282",
                                        "rs3845811", "rs4834792", "rs7041363", "rs7796", "rs8079811", "rs961764", 
                                        "rs9848170", "rs7310615"))]
write.csv(snps, file = "./result/MVMR_ALT-HbA1c-SBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_ALT-HbA1c-SBP_RETINA/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:4)], YG[,2], XGs_se[,c(2:4)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_ALT-HbA1c-SBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_ALT-HbA1c-SBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_ALT-HbA1c-SBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_ALT-HbA1c-SBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_ALT-HbA1c-SBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_ALT-HbA1c-SBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
library(robustMVMR)
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:4)])
sebetaGX = as.matrix(XGs_se[,c(2:4)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:4)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=14.11, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=14.11, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-ALT-HbA1c-DBP/RETINA2 #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/MVMR_ALT-HbA1c-DBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30620_raw", "ebi-a-GCST90002244", "ieu-b-39")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_ALT-HbA1c-DBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_ALT-HbA1c-DBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs11664194", "rs12321", "rs13395911", "rs1528293", "rs2160236", 
                                        "rs28667801", "rs2954027", "rs3802517", "rs4782568", "rs61912333", 
                                        "rs7041363", "rs710249", "rs7694000", "rs7938342", "rs9889262", "rs9893005"))]
write.csv(snps, file = "./result/MVMR_ALT-HbA1c-DBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_ALT-HbA1c-DBP_RETINA/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:4)], YG[,2], XGs_se[,c(2:4)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_ALT-HbA1c-DBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_ALT-HbA1c-DBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_ALT-HbA1c-DBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_ALT-HbA1c-DBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_ALT-HbA1c-DBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_ALT-HbA1c-DBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
library(robustMVMR)
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:4)])
sebetaGX = as.matrix(XGs_se[,c(2:4)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:4)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=14.11, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=14.11, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-AST-HbA1c-SBP/RETINA2 #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/MVMR_AST-HbA1c-SBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30650_raw", "ebi-a-GCST90002244", "ieu-b-38")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_AST-HbA1c-SBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_AST-HbA1c-SBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs1012089", "rs12321", "rs1821002", "rs1870735", "rs2024385",
                                        "rs3802517", "rs3828282", "rs3845811", "rs4834792", "rs7310615", "rs7463212", 
                                        "rs7796", "rs8030856", "rs8079811", "rs961764", "rs9848170"))]
write.csv(snps, file = "./result/MVMR_AST-HbA1c-SBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_AST-HbA1c-SBP_RETINA/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:4)], YG[,2], XGs_se[,c(2:4)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_AST-HbA1c-SBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_AST-HbA1c-SBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_AST-HbA1c-SBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_AST-HbA1c-SBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_AST-HbA1c-SBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_AST-HbA1c-SBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
library(robustMVMR)
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:4)])
sebetaGX = as.matrix(XGs_se[,c(2:4)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:4)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=10.65, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=10.65, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-AST-HbA1c-DBP/RETINA2 #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
rm(list = ls())
gc()
dir.create("./result/MVMR_AST-HbA1c-DBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30650_raw", "ebi-a-GCST90002244", "ieu-b-39")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_AST-HbA1c-DBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_AST-HbA1c-DBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs11664194", "rs12321", "rs1528293", "rs17696749", "rs2160236", 
                                        "rs28667801", "rs3802517", "rs4909314", "rs61912333", "rs710249", 
                                        "rs7694000", "rs7938342", "rs9889262", "rs9893005", "rs990619"))]
write.csv(snps, file = "./result/MVMR_AST-HbA1c-DBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_AST-HbA1c-DBP_RETINA/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:4)], YG[,2], XGs_se[,c(2:4)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_AST-HbA1c-DBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_AST-HbA1c-DBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_AST-HbA1c-DBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_AST-HbA1c-DBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_AST-HbA1c-DBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3], XGs_betas[,4]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3], XGs_se[,4]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_AST-HbA1c-DBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
library(robustMVMR)
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:4)])
sebetaGX = as.matrix(XGs_se[,c(2:4)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:4)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_AST-HbA1c-DBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=10.65, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=10.65, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-ALT/RETINA2(ajusted for HbA1C) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_ALT-HbA1c_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30620_raw", "ebi-a-GCST90002244")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_ALT-HbA1c_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_ALT-HbA1c_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs13395911", "rs2954027", "rs4782568", "rs7041363"))]
write.csv(snps, file = "./result/MVMR_ALT-HbA1c_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_ALT-HbA1c_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_ALT-HbA1c_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_ALT-HbA1c_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_ALT-HbA1c_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_ALT-HbA1c_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-HbA1c_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=14.11, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=14.11, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-ALT/RETINA2(ajusted for SBP) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_ALT-SBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30620_raw", "ieu-b-38")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_ALT-SBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_ALT-SBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs1012089", "rs12321", "rs13395911", "rs1821002", "rs1870735", 
                                        "rs2024385", "rs2954027", "rs3802517", "rs3828282", "rs3845811",
                                        "rs4834792", "rs7041363", "rs7310615", "rs7796", "rs8079811", 
                                        "rs961764", "rs9848170"))]
write.csv(snps, file = "./result/MVMR_ALT-SBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_ALT-SBP_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_ALT-SBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_ALT-SBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_ALT-SBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_ALT-SBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=14.11, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=14.11, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-ALT/RETINA2(ajusted for DBP) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_ALT-DBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30620_raw", "ieu-b-39")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_ALT-DBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_ALT-DBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs11664194", "rs12321", "rs13395911", "rs1528293", "rs17696749",
                                        "rs2160236", "rs28667801", "rs2954027", "rs3802517", "rs61912333",
                                        "rs7041363", "rs710249", "rs7694000", "rs7938342", "rs9889262", 
                                        "rs9893005"))]
write.csv(snps, file = "./result/MVMR_ALT-DBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_ALT-SBP_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_ALT-DBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_ALT-DBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_ALT-DBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_ALT-DBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_ALT-DBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_ALT-DBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=14.11, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=14.11, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-AST/RETINA2(ajusted for HbA1C) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_AST-HbA1C_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30650_raw", "ebi-a-GCST90002244")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_AST-HbA1C_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_AST-HbA1C_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs4924417", "rs13395911", "rs340882", "rs5112"))]
write.csv(snps, file = "./result/MVMR_AST-HbA1C_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_AST-HbA1C_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_AST-HbA1C_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_AST-HbA1C_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_AST-HbA1C_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_AST-HbA1c_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_AST-HbA1C_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_AST-HbA1C_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=10.65, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=10.65, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-AST/RETINA2(ajusted for SBP) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_AST-SBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30650_raw", "ieu-b-38")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_AST-SBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_AST-SBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs1012089", "rs12321", "rs1821002", "rs1870735", "rs2024385",
                                        "rs3802517", "rs3828282", "rs3845811", "rs4834792", "rs7310615",
                                        "rs7463212", "rs7796", "rs8030856", "rs8079811", "rs961764", 
                                        "rs9848170"))]
write.csv(snps, file = "./result/MVMR_AST-SBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_AST-SBP_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_AST-SBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_AST-SBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_AST-SBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_AST-SBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_AST-SBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_AST-SBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=10.65, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=10.65, p=mr_mvegger@Pvalue.Est[1])

############################### MVMR-AST/RETINA2(ajusted for DBP) #####################################
#加载R包
library(TwoSampleMR)
library(MRInstruments)
library(MVMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(tidyverse)
library(robustMVMR)
rm(list = ls())
gc()
dir.create("./result/MVMR_AST-DBP_RETINA2")

#读取exposure数据IVs
id_exposure = c("ukb-d-30650_raw", "ieu-b-39")
exposure = mv_extract_exposures(id_exposure = id_exposure)

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
snps = unique(exposure$SNP)
outcome_raw = outcome_raw[outcome_raw$SNP %in% snps,]
write.csv(outcome_raw, "./result/MVMR_AST-DBP_RETINA2/outcome.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/MVMR_AST-DBP_RETINA2/outcome.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
snps[!(snps %in% outcome$SNP)]
outcome$SNP[outcome$pval.outcome < 5e-08]
outcome = outcome[outcome$pval.outcome > 5e-08,]

#数据融合
data = mv_harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
snps = outcome$SNP[!(outcome$SNP %in% c("rs11664194", "rs12321", "rs1528293", "rs17696749", "rs2160236",
                                        "rs28667801", "rs3802517", "rs4909314", "rs61912333", "rs710249", 
                                        "rs7694000", "rs7938342", "rs9889262", "rs9893005"))]
write.csv(snps, file = "./result/MVMR_AST-DBP_RETINA2/snps.csv")

#去除与混淆因素相关的snp
excluded = read.table("./result/MVMR_AST-DBP_RETINA2/excluded.txt")
excluded = excluded$V1
snps = snps[!(snps %in% excluded)]
unique(excluded)

#数据整理
exposure = exposure[exposure$SNP %in% snps,]
outcome = outcome[outcome$SNP %in% snps,]

XGs_betas = exposure[,c("SNP", "beta.exposure", "exposure")]
XGs_betas = spread(XGs_betas, exposure, beta.exposure)
XGs_se = exposure[,c("SNP", "se.exposure", "exposure")]
XGs_se = spread(XGs_se, exposure, se.exposure)
YG = outcome[,c("SNP", "beta.outcome", "se.outcome")]

XGs_betas = XGs_betas[order(XGs_betas$SNP),]
XGs_se = XGs_se[order(XGs_se$SNP),]
YG = YG[order(YG$SNP),]
mvmr = format_mvmr(XGs_betas[,c(2:3)], YG[,2], XGs_se[,c(2:3)], YG[,3], XGs_betas$SNP)

#Instrument strength 
strength_mvmr = strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr = pleiotropy_mvmr(mvmr, gencov=0)
write.csv(strength_mvmr, file = "./result/MVMR_AST-DBP_RETINA2/F.csv")
save(pleiotropy_mvmr, file = "./result/MVMR_AST-DBP_RETINA2/pleiotropy.Rdata")

#MVMR-legacy
mvmr_res = mvmr(mvmr, 0, 1)
save(mvmr_res, file="./result/MVMR_AST-DBP_RETINA2/mvmr_res_final.Rdata")

#MVMR-IVW
mvmr_res_ivw = ivw_mvmr(mvmr)
save(mvmr_res_ivw, file="./result/MVMR_AST-DBP_RETINA2/mvmr_res_ivw_final.Rdata")

#MVMR-IVW
mr_mvivw = mr_mvivw(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                               bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]))
mr_mvivw
save(mr_mvivw, file="./result/MVMR_AST-DBP_RETINA2/mr_mvivw_final.Rdata")

#MVMR-Egger
mr_mvegger = mr_mvegger(mr_mvinput(bx = cbind(XGs_betas[,2], XGs_betas[,3]),
                                   bxse = cbind(XGs_se[,2], XGs_se[,3]), by = YG[,2], YG[,3]), orientate = 1)
mr_mvegger
save(mr_mvegger, file="./result/MVMR_AST-DBP_RETINA2/mr_mvegger_final.Rdata")

#Q-minimization
## -- SNP-outcome data
betaGY =  YG$beta.outcome
sebetaGY = YG$se.outcome
outcome = outcome[order(outcome$SNP),]
pvalbetaGY = outcome$pval.outcome
## -- SNP-exposure data
betaGX = as.matrix(XGs_betas[,c(2:3)])
sebetaGX = as.matrix(XGs_se[,c(2:3)])
XGs_p = exposure[,c("SNP", "pval.exposure", "exposure")]
XGs_p = spread(XGs_p, exposure, pval.exposure)
pvalbetaGX = as.matrix(XGs_p[,c(2:3)])
## -- Robust MVMR
fit = robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                 betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                 pval_threshold = 5e-08, plot = FALSE)
fit$rho_Exposures
mvmrcovmatrix = fit$rho_Exposures
mvmr_res_qhet = qhet_mvmr(mvmr, mvmrcovmatrix, T, 10000)
mvmr_res_qhet
save(mvmr_res_qhet, file="./result/MVMR_ALT-SBP_RETINA2/mvmr_res_qhet_final.Rdata")

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=mr_mvivw@Estimate[["Bx1"]], sd=10.65, p=mr_mvivw@Pvalue[["Bx1"]])
calculate_OR(beta=mr_mvegger@Estimate[1], sd=10.65, p=mr_mvegger@Pvalue.Est[1])

############################### Mediator effect ################################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(BWMR)
rm(list = ls())
gc()
dir.create("./result/Mediator")

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ebi-a-GCST90002244", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/Mediator/exposure_raw.csv")

#读取outcome数据
outcome = extract_outcome_data(snps = exposure$SNP, outcomes = "ieu-b-39", proxies = F)
write.csv(outcome, "./result/Mediator/outcome.csv")

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/Mediator/data_harmonised.csv")

#去除与混淆因素相关的SNP
data = data[!(data$SNP %in% snp_excluded),]

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/Mediator/MR_result.csv")
pdf(file = "./result/Mediator/scatterPlot1.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/Mediator/MR_heterogeneity1.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/Mediator/MR_pleiotropy1.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/Mediator/funnelPlot1.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/Mediator/leaveoneoutPlot1.pdf", width = 8, height = 16)
p2
dev.off()

#读取exposure数据IVs
exposure = extract_instruments(outcomes = "ieu-b-39", p1 = 5e-08, clump = TRUE,
                               r2 = 0.001, kb = 10000, access_token = NULL)
write.csv(exposure, file = "./result/Mediator/exposure_raw2.csv")

#读取outcome数据
outcome_raw = fread("E:/MR/finn_DM/finngen_R9_E4_DM2OPTH")
colnames(outcome_raw)
head(outcome_raw)
outcome_raw = outcome_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(outcome_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
outcome_raw = outcome_raw[outcome_raw$SNP %in% exposure$SNP,]
write.csv(outcome_raw, "./result/Mediator/outcome2.csv")

#将暴露SNP从结局中提取出来
merge_data = merge(exposure, outcome_raw, by = "SNP")

#读取outcome数据
outcome = read_outcome_data(snps = merge_data$SNP,
                            filename = "./result/Mediator/outcome2.csv",sep = ",",
                            snp_col = "SNP",beta_col = "beta",se_col = "se",
                            effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
exposure$SNP[!(exposure$SNP %in% outcome_raw$SNP)]

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]
write.csv(data, file = "./result/Mediator/data_harmonised2.csv")

#去除与混淆因素相关的SNP
data = data[!(data$SNP %in% snp_excluded),]

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result
write.csv(MR_result, file = "./result/Mediator/MR_result2.csv")
pdf(file = "./result/Mediator/scatterPlot2.pdf", width = 8, height = 8)
mr_scatter_plot(mr_results = mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),data)
dev.off()

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity
write.csv(mr_heterogeneity, file = "./result/Mediator/MR_heterogeneity2.csv")

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy
write.csv(mr_pleiotropy, file = "./result/Mediator/MR_pleiotropy2.csv")

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pdf(file = "./result/Mediator/funnelPlot2.pdf", width = 8, height = 8)
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
pdf(file = "./result/Mediator/leaveoneoutPlot2.pdf", width = 8, height = 40)
p2
dev.off()

# mediator-outcome 0.01316645(8.768948e-05 - 0.02624521)
total = 0.0289294120695198
exposure_mediator = 0.090083624
mediator_outcome = 0.01316645
mediation_effect = exposure_mediator * mediator_outcome
proportion = mediation_effect/total


############################### Bidirectional MR #########################
#打开R包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(BWMR)
rm(list = ls())
gc()
dir.create("./result/BUVMR_OPTH_ALT")

# 读取exposure数据IVs
exposure_raw = fread("E:/MR/finn_DM/finngen_R9_DM_RETINOPATHY_EXMORE")
exposure_raw = exposure_raw[,c("rsids", "beta", "sebeta", "af_alt", "alt", "ref", "pval")]
colnames(exposure_raw) = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p")
exposure_raw = exposure_raw[exposure_raw$p < 5e-08,]
exposure_raw = na.omit(exposure_raw)
exposure_raw = exposure_raw[!(str_detect(exposure_raw$SNP, ",")),]
exposure = format_data(exposure_raw, snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
                       effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p")
exposure = clump_data(exposure, clump_r2 = 0.001, clump_p1 = 5e-08, clump_p2 = 5e-08)

# 获取outcome数据
outcome = extract_outcome_data(snps = exposure$SNP, outcomes = "ukb-d-30650_raw", proxies = F)

#数据融合
data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
data$SNP[data$pval.outcome < 5e-08]
data = data[data$pval.outcome > 5e-08,]
data = data[data$mr_keep == TRUE,]

#MR PRESSO
mr_presso = run_mr_presso(data, NbDistribution = 1000)
mr_presso

#mr_method_list()
MR_result = generate_odds_ratios(mr_res = mr(data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
MR_result

#异质性检测
mr_heterogeneity = mr_heterogeneity(data)
mr_heterogeneity

#多效性检测
mr_pleiotropy = mr_pleiotropy_test(data)
mr_pleiotropy

#漏斗图
p1 = mr_funnel_plot(singlesnp_results = mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
p1
dev.off()

#留一法
p2 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
p2
dev.off()

#BWMR
library(BWMR)
BWMR_res = BWMR(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
BWMR_res

#OR计算(per SD)
calculate_OR = function(beta=NULL, sd=NULL, p=NULL){
  OR = exp(beta * sd)
  se = abs(log(OR)/qnorm(p/2))
  upper = OR + se * 1.96
  lower = OR - se * 1.96
  return(paste0(round(OR,3), "(", round(lower,3), "-", round(upper,3), ")"))
}
calculate_OR(beta=MR_result[1,7], sd=14.11, p=MR_result[1,9])


