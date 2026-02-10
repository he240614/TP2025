
library(tidyverse)
library(data.table)

setwd("D:/R/MDD_Metabolism")

rm(list = ls())
gc()

# rsid位置信息
pos = fread("snp_pos.txt.gz") %>% data.frame(check.names = F)
head(pos)

eqtl = fread("sc-eQTL/Bryois2022NN_OPCs_qtl.tsv") %>% data.frame(check.names = F)
head(eqtl)
pos1 = pos[pos$SNP %in% eqtl$variantId, c(1,4:6)]
# 合并rsid位置信息
eqtl = merge(pos1, eqtl, by.x = "SNP", by.y = "variantId")
# p值筛选
eqtl_pval = subset(eqtl, pValue < 5e-06)
gene = unique(eqtl_pval$geneName)

# 添加基因位置信息
grch37 = read.csv("GRCh37.csv", check.names = F)
grch37 = grch37[,c(1,2,12)]
anno = grch37[grch37$gene_name %in% gene,]
anno = anno[!duplicated(anno$gene_name),]
gene1 = gene[!gene %in% anno$gene_name]

grch38 = read.csv("GRCh38.csv", check.names = F)
grch38 = grch38[,c(1,2,12)]
anno1 = grch38[grch38$gene_name %in% gene1,]
anno = rbind(anno, anno1)

df = merge(eqtl_pval, anno, by.x = "geneName", by.y = "gene_name", all = T)
df$eaf = 1-df$MAF
df$probe = df$geneName
df$Orientation = "N"
head(df)
df1 = df %>% 
  dplyr::select(SNP,chrom,position,effect_allele,other_allele,eaf,
                probe,seqnames,start,geneName,Orientation,beta,se,pValue)
colnames(df1) = c("SNP","Chr","BP","A1","A2","Freq",
                  "Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","se","p")
df2 = df1 %>% filter(!is.na(SNP)) %>% 
  filter(!grepl(",", SNP)) %>% 
  filter(A1 %in% c("A","T","C","G")) %>% 
  filter(A2 %in% c("A","T","C","G"))
write.csv(df2, "sc-eQTL/OPCs_qtl.csv", row.names = F)

df2 = read.csv("sc-eQTL/OPCs_qtl.csv", check.names = F)
length(unique(df2$Probe))
write.table(df2, "sc-eQTL/OPCs_qtl.txt", sep = "\t", row.names = F, quote = F)

# txt转besd
# smr-1.3.1-win.exe --qfile ../OPCs_qtl.txt --make-besd --out ../OPCs/OPCs_qtl

# 结局数据处理
outcome = fread("D:/R/MDD_Metabolism/outcome/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz")
outcome = data.frame(outcome, check.names = F)
head(outcome)
outcome$samplesize = outcome$NCAS + outcome$NCON
outcome$eaf = (outcome$FCAS*outcome$NCAS+outcome$FCON*outcome$NCON)/outcome$samplesize
outcome$maf = ifelse(outcome$eaf<0.5, outcome$eaf, 1-outcome$eaf)

outcome_data = na.omit(outcome)
outcome_data$n = NA # 必须要有样本量这一列数据
head(outcome_data)
# 选取列(SNP名称,效应等位基因,非效应等位基因,效应等位基因频率,beta值,se值,p值,样本量)
outcome_data = dplyr::select(outcome_data, ID, EA, NEA, eaf, BETA, SE, PVAL, samplesize)
# 修改列名
colnames(outcome_data) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
head(outcome_data)
# 剔除SNP列为空白的行
outcome_data = subset(outcome_data, SNP != "")
# 剔除SNP列含有多个SNP的行
outcome_data = outcome_data[!grepl(",", outcome_data$SNP),]
# 去除重复的SNP
outcome_data = outcome_data %>% distinct(SNP, .keep_all = T)
# 保存数据
write.table(outcome_data, "D:/R/MDD_Metabolism/outcome/pgc_mdd.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# SMR分析
# smr-1.3.1-win.exe --bfile ../g1000_eur/g1000_eur --gwas-summary ../pgc_mdd.txt --beqtl-summary ../OPCs/OPCs_qtl --peqtl-smr 5e-6 --out OPCs_mdd --thread-num 10 --diff-freq-prop 0.9


# 读取结果数据
smr_Astros = read.delim("smr/result/Astros_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Astros$celltype = "Astros"
smr_Endo = read.delim("smr/result/Endo_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Endo$celltype = "Endo"
smr_Ex = read.delim("smr/result/Ex_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Ex$celltype = "Ex"
smr_Inhib = read.delim("smr/result/Inhib_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Inhib$celltype = "Inhib"
smr_Micro = read.delim("smr/result/Micro_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Micro$celltype = "Micro"
smr_Oligos = read.delim("smr/result/Oligos_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_Oligos$celltype = "Oligos"
smr_OPCs = read.delim("smr/result/OPCs_mdd.smr", header = TRUE, stringsAsFactors = FALSE)
smr_OPCs$celltype = "OPCs"

smr = rbind(smr_Astros,smr_Endo,smr_Ex,smr_Inhib,smr_Micro,smr_Oligos,smr_OPCs)
table(smr$celltype)

smr = smr[order(smr$p_SMR),]
smr$FDR = p.adjust(smr$p_SMR, method = "BH")
write.csv(smr, "smr/result/smr_all_fdr.csv", row.names = F)

