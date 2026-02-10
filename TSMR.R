
library(R.utils)
library(readr)
library(tidyr)
library(data.table)
library(dplyr)
library(vroom)
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
library(ieugwasr)
library(plinkbinr)
library(TwoSampleMR)
library(tidyverse)
library(future)
library(furrr)
library(ggvenn)
library(VennDiagram)
library(forestploter)
library(coloc)
library(openxlsx)
library(locuscomparer)

rm(list = ls())
gc()

# 结局数据处理 -----
setwd("D:/R/MDD_Metabolism")
outcome = fread("outcome/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz")
outcome = data.frame(outcome,check.names = F)
head(outcome)
outcome$samplesize = outcome$NCAS + outcome$NCON
outcome$eaf = (outcome$FCAS*outcome$NCAS+outcome$FCON*outcome$NCON)/outcome$samplesize
outcome$maf = ifelse(outcome$eaf<0.5, outcome$eaf, 1-outcome$eaf)
outcome = outcome[,c("#CHROM","POS","ID","EA","NEA","BETA","SE","PVAL","samplesize","NCAS","NCON","eaf","maf")]
saveRDS(outcome, file = "outcome/PGC_MDD.rds")

# 循环MR分析 -----
# 结局数据读取
outcome = readRDS("outcome/PGC_MDD.rds")
# rsid位置信息
pos = fread("snp_pos.txt.gz") %>% data.frame(check.names = F)
head(pos)

# 暴露数据读取
eqtl = read_tsv("D:/R/MDD_Metabolism/sc-eQTL/Bryois2022NN_Astros_qtl.tsv") %>% data.frame(check.names = F)
pos1 = pos[pos$SNP %in% eqtl$variantId, c(1,4:6)]
# 合并rsid位置信息
eqtl = merge(pos1, eqtl, by.x = "SNP", by.y = "variantId")
# 添加样本量
eqtl$samplesize = 2750
# p值筛选
eqtl_pval = subset(eqtl, pValue < 5e-06)
# 提取基因名
list_Probe = unique(eqtl_pval$geneName)

# MR分析
options(future.globals.maxSize = 10*1024^3)
plan(multicore, workers = 6)
list_Probe %>% future_map(~{
  exposure = subset(eqtl_pval, geneName == .x)
  exposure$eaf = 1-exposure$MAF
  tryCatch({
    exposure = format_data(exposure,
                           type = "exposure",
                           chr_col = "chrom",
                           pos_col = "position",
                           snp_col = "SNP",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           beta_col = "beta",
                           se_col = "se",
                           pval_col = "pValue",
                           eaf_col = "eaf",
                           samplesize_col = "samplesize",
                           id_col = "geneName",
                           phenotype_col = "cellTypeName")
    # 本地去除连锁不平衡
    exposure_clump = ieugwasr::ld_clump(dplyr::tibble(rsid = exposure$SNP,
                                                      pval = exposure$pval.exposure),
                                        clump_kb = 10000,
                                        clump_r2 = 0.001,
                                        clump_p = 1,
                                        bfile = "D:/R/MR-Class/TSMR/1kg.v3/EUR",
                                        plink_bin = plinkbinr::get_plink_exe(),
                                        pop = "EUR")
    exposure_clump = exposure[exposure$SNP %in% exposure_clump$rsid,]
    # 计算F值
    R2a = 2*(exposure_clump$beta.exposure^2)*exposure_clump$eaf.exposure*(1-exposure_clump$eaf.exposure)
    R2b = 2*(exposure_clump$se.exposure^2)*exposure_clump$samplesize.exposure*exposure_clump$eaf.exposure*(1-exposure_clump$eaf.exposure)
    R2 = R2a/(R2a+R2b)
    exposure_clump$Fz = R2*(exposure_clump$samplesize.exposure-2)/(1-R2)
    # 结局数据处理
    outcome_dat = outcome[outcome$ID %in% exposure_clump$SNP,]
    outcome_dat = format_data(dat = outcome_dat,
                              type = "outcome",
                              chr_col = "#CHROM",
                              pos_col = "POS",
                              snp_col = "ID",
                              effect_allele_col = "EA",
                              other_allele_col = "NEA",
                              beta_col = "BETA",
                              pval_col = "PVAL",
                              se_col = "SE",
                              eaf_col = "eaf",
                              samplesize_col = "samplesize")
    outcome_dat$id.outcome = "MDD"
    outcome_dat = subset(outcome_dat, pval.outcome>5e-06)
    # 数据整合
    mr_data = harmonise_data(exposure_dat = exposure_clump,
                             outcome_dat = outcome_dat,
                             action = 2)
    write.csv(mr_data, paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_data/",.x,"_data.csv"))
    # MR分析
    res = mr(mr_data)
    res = generate_odds_ratios(res)
    write.csv(res, paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_result/",.x,"_res.csv"))
    # 敏感性分析
    het = mr_heterogeneity(mr_data)
    write.csv(het, paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_heterogeneity/",.x,"_het.csv"))
    pleio = mr_pleiotropy_test(mr_data)
    write.csv(pleio, paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_pleiotropy/",.x,"_pleio.csv"))
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
})

result = NULL
heterogeneity = NULL
pleiotropy = NULL
for (i in 1:length(list_Probe)) {
  tryCatch({
    res = read.csv(paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_result/",list_Probe[i],"_res.csv"), row.names = 1, check.names = F)
    result = rbind(result, res)
    het = read.csv(paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_heterogeneity/",list_Probe[i],"_het.csv"), row.names = 1, check.names = F)
    heterogeneity = rbind(heterogeneity, het)
    pleio = read.csv(paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_pleiotropy/",list_Probe[i],"_pleio.csv"), row.names = 1, check.names = F)
    pleiotropy = rbind(pleiotropy, pleio)
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}
write.csv(result, paste0("D:/R/MDD_Metabolism/PGC/Astros/MR_result.csv"))
write.csv(heterogeneity, "D:/R/MDD_Metabolism/PGC/Astros/MR_heterogeneity.csv")
write.csv(pleiotropy, "D:/R/MDD_Metabolism/PGC/Astros/MR_pleiotropy.csv")

# nohup Rscript MR.R &

# mr_data合并
# 暴露数据读取
eqtl = read_tsv("D:/R/MDD_Metabolism/sc-eQTL/Bryois2022NN_Oligos_qtl.tsv") %>% data.frame(check.names = F)
pos1 = pos[pos$SNP %in% eqtl$variantId, c(1,4:6)]
# 合并rsid位置信息
eqtl = merge(pos1, eqtl, by.x = "SNP", by.y = "variantId")
# 添加样本量
eqtl$samplesize = 2750
# p值筛选
eqtl_pval = subset(eqtl, pValue < 5e-06)
# 提取基因名
list_Probe = unique(eqtl_pval$geneName)

data = NULL
for (i in 1:length(list_Probe)) {
  tryCatch({
    mr_data = read.csv(paste0("D:/R/MDD_Metabolism/PGC/Oligos/MR_data/",list_Probe[i],"_data.csv"),
                       row.names = 1, check.names = F)
    data = rbind(data, mr_data)
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}
data$effect_allele.exposure = gsub("TRUE", "T", data$effect_allele.exposure)
data$other_allele.exposure = gsub("TRUE", "T", data$other_allele.exposure)
data$effect_allele.outcome = gsub("TRUE", "T", data$effect_allele.outcome)
data$other_allele.outcome = gsub("TRUE", "T", data$other_allele.outcome)
write.csv(data, "D:/R/MDD_Metabolism/PGC/Oligos/MR_data.csv", row.names = F)


# 结果筛选 -----
cell = c("Astros","Endo","Ex","Inhib","Micro","Oligos","OPCs")
for (i in 1:length(cell)) {
  # i = 1
  print(i)
  tryCatch({
    result = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/MR_result.csv"), row.names = 1, check.names = F)
    # 提取IVW和Wald ratio分析结果
    res = result[result$method %in% c("Wald ratio","Inverse variance weighted"),]
    # 计算FDR
    res = res[order(res$pval),]
    res$fdr = p.adjust(res$pval, method = "BH")
    write.csv(res, paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/MR_result_fdr.csv"))
    # 筛选阳性结果
    res_fdr = res[res$fdr < 0.05,]
    res_p = res[res$pval < 0.05,]
    # 结果保存
    write.csv(res_fdr, paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/MR_result_fdr_p0.05.csv"))
    write.csv(res_p, paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/MR_result_p0.05.csv"))
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}

cell = c("Astros","Endo","Ex","Inhib","Micro","Oligos","OPCs")
all_res = NULL
for (i in 1:length(cell)) {
  print(i)
  # i = 1
  res = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/MR_result.csv"), row.names = 1, check.names = F)
  res$celltype = cell[i]
  all_res = rbind(all_res, res)
}
table(all_res$celltype)
all_res = all_res[order(all_res$pval),]
all_res$FDR = p.adjust(all_res$pval, method = "BH")
write.csv(all_res, "PGC/all_MRresult_fdr.csv", row.names = F)



# Steiger方向检验 -----
cell = c("Astros","Endo","Ex","Inhib","Micro","Oligos","OPCs")
j = 7
res = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[j],"/MR_result.csv"), row.names = 1, check.names = F)
symbol = res$id.exposure
Steiger = NULL
for (i in 1:length(symbol)) {
  print(i)
  # i = 1
  mr_data = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[j],"/MR_data/",symbol[i],"_data.csv"), row.names = 1, check.names = F)
  out = directionality_test(mr_data)
  Steiger = rbind(Steiger, out)
}
write.csv(Steiger, paste0("D:/R/MDD_Metabolism/PGC/",cell[j],"/MR_Steiger.csv"))


# MR结果与差异分析结果取交集 -----
setwd("D:/R/MDD_Metabolism/wilcox")
sce_diff = read.csv("wilcox.allcell.csv", row.names = 1, check.names = F)
cell = unique(sce_diff$cell)
for (i in 1:length(cell)) {
  print(i)
  # i = 1
  cell[i]
  df = sce_diff[sce_diff$cell == cell[i],]
  df = df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.2,]
  # df = df[df$p_val_adj < 0.05,]
  df_up = df[df$avg_log2FC > 0,]
  df_down = df[df$avg_log2FC < 0,]
  mr_res = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/Pval/MR_result_p0.05.csv"), row.names = 1, check.names = F)
  mr_risk = mr_res[mr_res$or > 1,]
  mr_pro = mr_res[mr_res$or < 1,]
  
  venn_list_up = list(diff_up = df_up$gene_name, MR_risk = mr_risk$id.exposure)
  names(venn_list_up) = c(paste0(cell[i],"_diff_up"), paste0(cell[i],"_MR_risk"))
  p1 = ggvenn(data = venn_list_up,         #数据列表
              columns = NULL,            #对选中的列名绘图，最多选择4个，NULL为默认全选
              show_elements = F,         #当为TRUE时，显示具体的交集情况，而不是交集个数
              label_sep = "\n",          #当show_elements = T时生效，分隔符 \n 表示的是回车的意思
              show_percentage = F,       #显示每一组的百分比
              digits = 1,                #百分比的小数点位数
              fill_color = c("#db6968","#4d97cd"), #填充颜色
              fill_alpha = 0.6,          #填充透明度
              stroke_color = "black",    #边缘颜色
              stroke_alpha = 0.8,        #边缘透明度
              stroke_size = 1,           #边缘粗细
              stroke_linetype = "twodash", #边缘线条  实线：solid  虚线：twodash longdash  点：dotdash dotted dashed  无：blank
              set_name_color = c("#db6968","#4d97cd"),  #组名颜色
              set_name_size = 8,         #组名大小
              text_color = "black",      #交集个数颜色
              text_size = 9)+            #交集个数文字大小
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  p1
  ggsave(paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_Venn_up.pdf"), width = 8, height = 8, plot = p1)
  ggsave(paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_Venn_up.tiff"), width = 8, height = 8, plot = p1, dpi = 300)
  
  venn_list_down = list(diff_down = df_down$gene_name, MR_pro = mr_pro$id.exposure)
  names(venn_list_down) = c(paste0(cell[i],"_diff_down"), paste0(cell[i],"_MR_pro"))
  p2 = ggvenn(data = venn_list_down,         #数据列表
              columns = NULL,            #对选中的列名绘图，最多选择4个，NULL为默认全选
              show_elements = F,         #当为TRUE时，显示具体的交集情况，而不是交集个数
              label_sep = "\n",          #当show_elements = T时生效，分隔符 \n 表示的是回车的意思
              show_percentage = F,       #显示每一组的百分比
              digits = 1,                #百分比的小数点位数
              fill_color = c("#459943","#f8984e"), #填充颜色
              fill_alpha = 0.6,          #填充透明度
              stroke_color = "black",    #边缘颜色
              stroke_alpha = 0.8,          #边缘透明度
              stroke_size = 1,           #边缘粗细
              stroke_linetype = "twodash", #边缘线条  实线：solid  虚线：twodash longdash  点：dotdash dotted dashed  无：blank
              set_name_color = c("#459943","#f8984e"),  #组名颜色
              set_name_size = 8,         #组名大小
              text_color = "black",      #交集个数颜色
              text_size = 9)+             #交集个数文字大小
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  p2
  ggsave(paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_Venn_down.pdf"), width = 8, height = 8, plot = p2)
  ggsave(paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_Venn_down.tiff"), width = 8, height = 8, plot = p2, dpi = 300)
  # 提取交集基因
  inter_up = get.venn.partitions(venn_list_up)
  gene_up = as.data.frame(inter_up[1,4]) %>% mutate(condition = "Up")
  inter_down = get.venn.partitions(venn_list_down)
  gene_down = as.data.frame(inter_down[1,4]) %>% mutate(condition = "Down")
  gene = rbind(gene_up, gene_down)
  write.csv(gene, paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_inter_gene_wilcox.csv"), row.names = F)
}

# 森林图 -----
setwd("D:/R/MDD_Metabolism/wilcox")
# 交集基因结果合并
sce_diff = read.csv("wilcox.allcell.csv", row.names = 1, check.names = F)
cell = unique(sce_diff$cell)
rb = NULL
for (i in 1:length(cell)) {
  print(i)
  # i = 1
  cell[i]
  df = sce_diff[sce_diff$cell == cell[i],]
  tryCatch({
    mr_res = read.csv(paste0("D:/R/MDD_Metabolism/PGC/",cell[i],"/Pval/MR_result_p0.05.csv"), row.names = 1, check.names = F)
    gene = read.csv(paste0("D:/R/MDD_Metabolism/logFC_0.2/",cell[i],"_inter_gene_wilcox.csv"), check.names = F)
    
    MR = mr_res[mr_res$id.exposure %in% gene$X1,] %>% 
      dplyr::select(id.exposure,id.outcome,method,nsnp,or,or_lci95,or_uci95,pval) %>% 
      dplyr::rename(Symbol = id.exposure, Outcome = id.outcome, Method = method, SNPs = nsnp, MR.Pval = pval)
    diff = df[df$gene_name %in% gene$X1,] %>% 
      dplyr::select(gene_name,avg_log2FC,p_val_adj,cell) %>% 
      dplyr::rename(Symbol = gene_name, log2FC = avg_log2FC, scRNA.Pval = p_val_adj, celltype = cell)
    merge = merge(MR, diff, by = "Symbol")
    rb = rbind(rb,merge)
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}
rb = rb %>% arrange(celltype, or)
rb$Method[rb$Method == "Inverse variance weighted"] = "IVW"
write.csv(rb, "D:/R/MDD_Metabolism/logFC_0.2/MR_diff_wilcox.csv", row.names = F)

OR = read.csv("MR_diff_venn.csv", check.names = F)
OR$MR.Pval = ifelse(OR$MR.Pval < 0.001, "<0.001", sprintf("%.3f", OR$MR.Pval))
OR$scRNA.Pval = ifelse(OR$scRNA.Pval < 0.001, "<0.001", sprintf("%.3f", OR$scRNA.Pval))
OR$log2FC = ifelse(is.na(OR$log2FC), "", sprintf("%.3f", OR$log2FC))
OR$`OR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", OR$or, OR$or_lci95, OR$or_uci95)
OR$` ` = paste(rep(" ", 20), collapse = " ")
OR$SNPs[is.na(OR$SNPs)] = " "
OR$MR.Pval[is.na(OR$MR.Pval)] = " "
OR$scRNA.Pval[is.na(OR$scRNA.Pval)] = " "
OR$`OR (95% CI)`[grepl("NA", OR$`OR (95% CI)`)] = " "

# OR$Symbol[c(2:5,7,9:21,23:24,26:43,45)] = paste0("     ", OR$Symbol[c(2:5,7,9:21,23:24,26:43,45)])
# OR$Symbol[c(2:7,9,11:22,24,26:43,45:46)] = paste0("     ", OR$Symbol[c(2:7,9,11:22,24,26:43,45:46)])
# OR$Symbol[c(2:12,14:16,18:29,31:34,36:57,59:63)] = paste0("     ", OR$Symbol[c(2:12,14:16,18:29,31:34,36:57,59:63)])
OR$Symbol[c(2:14,16:18,20:31,33:36,38:58,60:64)] = paste0("     ", OR$Symbol[c(2:14,16:18,20:31,33:36,38:58,60:64)])
OR$Outcome = paste0("    ", OR$Outcome)
OR$Method[which(nchar(OR$Method) == 3)] = paste0("    ", OR$Method[which(nchar(OR$Method) == 3)])
OR$SNPs[which(nchar(OR$SNPs) == 1)] = paste0("   ", OR$SNPs[which(nchar(OR$SNPs) == 1)])
OR$MR.Pval[which(nchar(OR$MR.Pval) == 5)] = paste0("   ", OR$MR.Pval[which(nchar(OR$MR.Pval) == 5)])
OR$MR.Pval[which(nchar(OR$MR.Pval) == 6)] = paste0(" ", OR$MR.Pval[which(nchar(OR$MR.Pval) == 6)])
OR$scRNA.Pval[which(nchar(OR$scRNA.Pval) == 5)] = paste0("       ", OR$scRNA.Pval[which(nchar(OR$scRNA.Pval) == 5)])
OR$scRNA.Pval[which(nchar(OR$scRNA.Pval) == 6)] = paste0("     ", OR$scRNA.Pval[which(nchar(OR$scRNA.Pval) == 6)])
OR$log2FC[which(nchar(OR$log2FC) == 5)] = paste0("  ", OR$log2FC[which(nchar(OR$log2FC) == 5)])
OR$log2FC[which(nchar(OR$log2FC) == 6)] = paste0(" ", OR$log2FC[which(nchar(OR$log2FC) == 6)])
colnames(OR)[1] = paste0("     ", colnames(OR)[1])
colnames(OR)[3] = paste0(" ", colnames(OR)[3])
colnames(OR)[11] = paste0("       ", colnames(OR)[11])

tm = forest_theme(base_size = 10,
                  title_gp = gpar(cex = 1.2, fontface = "bold"), # 文本的大小
                  ci_pch = 23,      # 可信区间点的形状
                  ci_col = "black", # CI的颜色
                  ci_fill = "red",  # CI中se点的颜色填充
                  ci_alpha = 0.8,   # CI透明度
                  ci_lty = 1,       # CI的线型
                  ci_lwd = 2,       # CI的线宽
                  ci_Theight = 0.2, # CI的高度，默认是NULL
                  refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey20"), #参考线
                  core = list(bg_params = list(fill = c("white")))) #设置背景为白色
p1 = forest(OR[,c(1:4,12,8,11,9,10)],
            est = OR$or, #效应值
            lower = OR$or_lci95, # 置信区间下限
            upper = OR$or_uci95, # 置信区间上限
            ci_column = 5, # 在哪一列画森林图，选空的那一列
            ref_line = 1, # 参考线位置
            xlim = c(0.92,1.08), # 设置轴范围
            ticks_at = c(0.94,0.97,1,1.03,1.06),# 设置刻度
            theme = tm)
p2 = edit_plot(p1, row = which(OR$or < 1), col = 5, which = "ci", gp = gpar(fill = "#105BA2", col = "grey40", cex = 1.2))
p3 = edit_plot(p2, row = which(OR$or > 1), col = 5, which = "ci", gp = gpar(fill = "#c6133b", col = "grey40", cex = 1.2))
p4 = add_border(p3, part = "header", where = "bottom", gp = gpar(lwd = 1.5))
p5 = add_border(p4, part = "header", where = "top", gp = gpar(lwd = 1.5))
p6 = add_border(p5, row = c(14,18,31,36,58,nrow(OR)), where = "bottom", gp = gpar(lwd = 1.5))
p7 = edit_plot(p6, row = c(1,15,19,32,37,59), gp = gpar(fontface = "bold", fontsize = 12))
# p6 = add_border(p5, row = c(12,16,29,34,57,nrow(OR)), where = "bottom", gp = gpar(lwd = 1.5))
# p7 = edit_plot(p6, row = c(1,13,17,30,35,58), gp = gpar(fontface = "bold", fontsize = 12))
# p6 = add_border(p5, row = c(5,7,21,24,43,nrow(OR)), where = "bottom", gp = gpar(lwd = 1.5))
# p7 = edit_plot(p6, row = c(1,6,8,22,25,44), gp = gpar(fontface = "bold", fontsize = 12))
# p6 = add_border(p5, row = c(7,9,22,24,43,nrow(OR)), where = "bottom", gp = gpar(lwd = 1.5))
# p7 = edit_plot(p6, row = c(1,8,10,23,25,44), gp = gpar(fontface = "bold", fontsize = 12))
p7
ggsave(paste0("graph_Pval/MR_forest.pdf"), width = 9, height = 15, plot = p7, bg = "white")
ggsave(paste0("graph_Pval/MR_forest.tiff"), width = 9, height = 15, plot = p7, dpi = 300, bg = "white")


# eqtl-MDD共定位分析 -----
setwd("D:/R/MDD_Metabolism")
# 读取整理MDD数据
outcome = readRDS("outcome/PGC_MDD.rds")
head(outcome)
gwas = outcome %>% dplyr::select("ID","#CHROM","POS","EA","NEA","maf","BETA","SE","PVAL","samplesize","NCAS")
colnames(gwas) = c("SNP","chrom","pos","effect_allele","other_allele","maf","beta","se","pval","samplesize","ncase")
# 计算共定位分析需要的数据
gwas$varbeta = gwas$se^2
gwas$s = gwas$ncase/gwas$samplesize
gwas$z = gwas$beta/gwas$se
gwas = subset(gwas, !duplicated(SNP))
gwas = na.omit(gwas)

# MR筛选结果
# mr_diff = read.xlsx("smr/list.xlsx", sheet = 1)
mr_diff = read.csv("smr/gene.csv", check.names = F)

# rsid位置信息
pos = fread("snp_pos.txt.gz") %>% data.frame(check.names = F)
head(pos)
celltype = unique(mr_diff$cell.x)
all_coloc = NULL

j = 1
symbol = mr_diff$probeID[mr_diff$cell.x == celltype[j]]
# 读取eqtl数据
all_eqtl = read_tsv(paste0("D:/R/MDD_Metabolism/02.TSMR分析/sc-eQTL/Bryois2022NN_",celltype[j],"_qtl.tsv")) %>% data.frame(check.names = F)
pos1 = pos[pos$SNP %in% all_eqtl$variantId, c(1,4:6)]
# 合并rsid位置信息
all_eqtl = merge(pos1, all_eqtl, by.x = "SNP", by.y = "variantId")
# 添加样本量
all_eqtl$samplesize = 2750

coloc = NULL
for (i in 1:length(symbol)) {
  print(i)
  i = 2
  eqtl = subset(all_eqtl, geneName == symbol[i])
  # 读取整理eQTL
  dat = eqtl %>% dplyr::select("SNP","chrom","position","effect_allele","other_allele","MAF","beta","se","pValue","samplesize","geneName")
  # 计算共定位分析需要的数据
  dat$varbeta = dat$se^2
  dat$z = dat$beta/dat$se
  
  # 挑选lead SNP
  lead = dat[order(dat$pValue),]
  leadchr = as.numeric(lead$chrom[1])
  leadpos = as.numeric(lead$position[1])
  
  QTLdata = dat[dat$chrom == leadchr & dat$position > leadpos-100000 & dat$position < leadpos+100000,]
  QTLdata = na.omit(QTLdata)
  QTLdata = subset(QTLdata, !duplicated(SNP))
  QTLdata = QTLdata[QTLdata$pValue > 0,]
  tryCatch({
    sameSNP = intersect(QTLdata$SNP, gwas$SNP)
    QTLdata = QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
    GWASdata = gwas[gwas$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
    
    # 共定位分析
    result = coloc.abf(dataset1 = list(pvalues = GWASdata$pval,
                                       snp = GWASdata$SNP,
                                       type = "cc",
                                       s = GWASdata$s[1],
                                       N = GWASdata$samplesize[1]),
                       dataset2 = list(pvalues = QTLdata$pValue,
                                       snp = QTLdata$SNP,
                                       type = "quant",
                                       N = QTLdata$samplesize[1]),
                       MAF = QTLdata$MAF)
    # 结果提取
    res = t(result$summary) %>% data.frame()
    res$symbol = symbol[i]
    res$celltype = celltype[j]
    coloc = rbind(coloc, res)
  },error = function(e) {
    #处理错误的代码，可以选择打印错误信息等
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}
all_coloc = rbind(all_coloc, coloc)
write.csv(all_coloc, "smr/eqtl_MDD_coloc.csv", row.names = F)


# 共定位可视化
gwas_fn = GWASdata[,c("SNP","pval")] %>% dplyr::rename(rsid = SNP, pval = pval)
eqtl_fn = QTLdata[,c("SNP","pValue")] %>% dplyr::rename(rsid = SNP, pval = pValue)
# 绘图
pdf(paste0("03.eQTL-MDD共定位分析/locuscompare_",symbol[i],".pdf"), width = 6, height = 5)
locuscompare(in_fn1 = gwas_fn, 
             in_fn2 = eqtl_fn, 
             title1 = "MDD", 
             title2 = symbol[i])
dev.off()
tiff(paste0("03.eQTL-MDD共定位分析/locuscompare_",symbol[i],".tiff"), width = 1800, height = 1500, res = 300)
locuscompare(in_fn1 = gwas_fn, 
             in_fn2 = eqtl_fn, 
             title1 = "MDD", 
             title2 = symbol[i])
dev.off()

