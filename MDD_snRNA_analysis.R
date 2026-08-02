## Script:
## 3.MDD_snRNA_analysis.R
## Purpose:
## MDD snRNA-seq analysis

rm(list = ls())
library(Seurat)
library(Matrix)

## Project directories
project_dir <- "."      # change based on your own path

data_dir <- file.path(project_dir, "01.rawdata")
QC_dir <- file.path(data_dir, "02.QC")
cellanno_dir <- file.path(data_dir, "03.CellAnno")
DEanalysis_dir <- file.path(data_dir, "04.DEanalysis")
edgeR_dir <- file.path(DEanalysis_dir, "edgeR")
MASTRE_dir <- file.path(DEanalysis_dir, "MASTRE")
sgGSEA_dir <- file.path(project_dir, "05.singlegeneGSEA")
COQ8A_dir <- file.path(project_dir, "06.COQ8AEx")
DAanalysis_dir <- file.path(project_dir, "07.DAanalysis")

dir_list <- c(data_dir, QC_dir, cellanno_dir, DEanalysis_dir,
              edgeR_dir, MASTRE_dir, sgGSEA_dir, COQ8A_dir,
              DAanalysis_dir)
invisible(
  lapply(dir_list, dir.create,
         recursive = TRUE,
         showWarnings = FALSE))
#######Data reading########
genes_file <- "GSE213982_combined_counts_matrix_genes_rows.csv.gz"
genes <- read.csv(gzfile(file.path(data_dir, genes_file)), header = TRUE, stringsAsFactors = FALSE)
gene_names <- genes[, 1]  
cells_file <- "GSE213982_combined_counts_matrix_cells_columns.csv.gz"
cells <- read.csv(gzfile(file.path(data_dir, cells_file)), header = TRUE, stringsAsFactors = FALSE)
cell_barcodes <- cells[, 1]  
matrix_file <- "GSE213982_combined_counts_matrix.mtx.gz"
count_matrix <- readMM(file = gzfile(file.path(data_dir, matrix_file)))
stopifnot(length(gene_names) == nrow(count_matrix),
          length(cell_barcodes) == ncol(count_matrix))
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- cell_barcodes
gse213.sc <- CreateSeuratObject(
  counts = count_matrix,
  project = "brain",
  assay = "RNA",
  min.cells = 3,
  min.features =200)
print(gse213.sc)
saveRDS(gse213.sc,file.path(data_dir, "gse213.sc.rds"))

##########QC#######
library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(clustree)
library(decontX)
library(patchwork)
library(harmony)
gse213.sc <- readRDS(file.path(data_dir, "gse213.sc.rds"))
rownames(gse213.sc@assays$RNA@layers$counts) = Features(gse213.sc)
colnames(gse213.sc@assays$RNA@layers$counts) = Cells(gse213.sc)
table(gse213.sc$orig.ident)
cell_names <- colnames(gse213.sc)
head(cell_names)
samples <- sapply(strsplit(cell_names, split = "\\."), function(x) x[1])
gse213.sc[["sample"]] <- samples
unique(gse213.sc@meta.data[, "sample"])
group <- read.csv(file.path(data_dir, "group.csv"),header = T)
samples_in_seurat <- gse213.sc$sample
groups_mapped <- group$group[match(samples_in_seurat, group$sample)]
gse213.sc[["group"]] <- groups_mapped
head(gse213.sc@meta.data[, c("sample", "group")])
table(gse213.sc$sample, gse213.sc$group)
# delete "M24_2"
gse213.sc <- subset(
  x = gse213.sc,
  subset = sample != "M24_2"
)
table(gse213.sc$sample)
gse213.sc[["percent.mt"]] = PercentageFeatureSet(gse213.sc, pattern = "^MT-")
p1 <- VlnPlot(object = gse213.sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
ggsave(file.path(QC_dir, "01_QCbefore.pdf"), width = 15, height = 6, plot = p1)
ggsave(file.path(QC_dir, "01_QCbefore.tiff"), width = 15, height = 6, dpi = 300, plot = p1)

plot1 = FeatureScatter(gse213.sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(gse213.sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 = plot1 + plot2
p2
ggsave(file.path(QC_dir, "02_FeatureScatter.pdf"), width = 12, height = 6, plot = p2)
ggsave(file.path(QC_dir, "02_FeatureScatter.tiff"), width = 12, height = 6, dpi = 300, plot = p2)

gse213.sc = subset(gse213.sc, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 100000 & percent.mt < 10)
table(gse213.sc$orig.ident)
p3 = VlnPlot(gse213.sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3
ggsave(file.path(QC_dir, "03_QCafter.pdf"), width = 15, height = 6, plot = p3)
ggsave(file.path(QC_dir, "03_QCafter.tiff"), width = 15, height = 6, dpi = 300, plot = p3)

decontX_results = decontX::decontX(gse213.sc@assays$RNA@layers$counts)
gse213.sc$Contamination = decontX_results$contamination
rt = gse213.sc@meta.data
rt = rt[which(rt$Contamination < 0.3),]
gse213.sc = gse213.sc[,row.names(rt)]

options(future.globals.maxSize= 8388608000)
gse213.sc <- SCTransform(
  gse213.sc,
  assay = "RNA",
  vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), 
  variable.features.n = 3000,  
  verbose = TRUE
)

head(VariableFeatures(gse213.sc), 20)
plot1 <- VlnPlot(gse213.sc, features = head(VariableFeatures(gse213.sc), 10), ncol = 5)
plot1 + RotatedAxis()
gse213.sc <- RunPCA(gse213.sc, assay = "SCT", features = VariableFeatures(gse213.sc))
p4 <- DimPlot(object = gse213.sc, reduction = "pca", pt.size = 0.1, group.by = "sample",raster=FALSE)
p4
ggsave(file.path(QC_dir, "04_PCA.pdf"), width = 8, height = 6, plot = p4)
ggsave(file.path(QC_dir, "04_PCA.tiff"), width = 8, height = 6, dpi = 300, plot = p4)
gse213.sc = RunHarmony(gse213.sc, group.by.vars = "sample", plot_convergence = TRUE)
p5 = DimPlot(object = gse213.sc, reduction = "harmony", pt.size = 0.1, group.by = "sample",raster=FALSE)
p5
ggsave(file.path(QC_dir, "05_Harmony.pdf"), width = 8, height = 6, plot = p5)
ggsave(file.path(QC_dir, "05_Harmony.tiff"), width = 8, height = 6, dpi = 300, plot = p5)
p6 = ElbowPlot(gse213.sc, ndims = 50, reduction = "pca")+theme_bw()
p6
ggsave(file.path(QC_dir, "06_Elbowplot.pdf"), width = 8, height = 6, plot = p6)
ggsave(file.path(QC_dir, "06_Elbowplot.tiff"), width = 8, height = 6, dpi = 300, plot = p6)
gse213.sc <- FindNeighbors(gse213.sc, assay = "SCT", reduction = "harmony",dims = 1:30)
gse213.sc = FindClusters(gse213.sc, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
p7 = clustree(gse213.sc)
p7
ggsave(file.path(QC_dir, "07_clustree.pdf"), width = 12, height = 10, plot = p7)
ggsave(file.path(QC_dir, "07_clustree.tiff"), width = 12, height = 10, dpi = 300, plot = p7)
Idents(gse213.sc) = "SCT_snn_res.0.6"
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead","#ffad60",
          "#7F95D1","#d9534f",  "#57D1C9", "#ED5485",  "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")
gse213.sc = RunUMAP(gse213.sc, reduction = "harmony", seed.use = 786452, dims = 1:20)
gse213.sc = RunTSNE(gse213.sc, reduction = "harmony", seed.use = 786452, dims = 1:20)
p8 = DimPlot(gse213.sc, reduction = "umap", group.by = "SCT_snn_res.0.6", cols = plan, label = T,repel = T,raster=FALSE)
p8
ggsave(file.path(QC_dir, "08_UMAP786452.pdf"), width = 8, height = 6, plot = p8)
ggsave(file.path(QC_dir, "08_UMAP786452.tiff"), width = 8, height = 6, dpi = 300, plot = p8)
save(gse213.sc, file = file.path(data_dir, "seurat_object_snn_0.6_1142.rdata"))

##########Cell annotation#########
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)
load(file.path(data_dir, "seurat_object_snn_0.6_1142.rdata"))
allmarkers = FindAllMarkers(gse213.sc, only.pos = TRUE, min.pct = 0.6, logfc.threshold = 0.5)
write.csv(allmarkers, file.path(cellanno_dir, "allmarkers.csv"))
Top10.coarse = allmarkers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(Top10.coarse, file.path(cellanno_dir, "Top10.markers.csv"))
####未设置随机种子
markers <- c(
  "NRGN","TUBA1B","SLC17A7",  #Excitatory Neurons(Ex) 
  "CERCAM","TF","ENPP2",   # Oligodendrocytes(Oligos) 
  "GAD1","GAD2", # Inhibitory Neurons (Inhib)  
  "ADGRV1","SLC1A2","FGFR3",  #Astrocytes (Astros) 
  "CLDN5","FLT1","ABCB1",  # Endothelial cells (Endo) 
  "P2RY12","CSF1R","DOCK8",   # Microglia (Micro) 
  "PDGFRA","VCAN","SEMA5A"  # OPCs (OPCs)  
)
p9 = DotPlot(gse213.sc, features = markers) + theme_bw() + RotatedAxis()
p9

gse213.sc$cell.types = recode(gse213.sc$SCT_snn_res.0.6,
                              "0" = "Oligos",
                              "1" = "Ex",
                              "2" = "Ex",
                              "3" = "Ex",
                              "4" = "Ex",
                              "5" = "Astros",
                              "6" = "Inhib",
                              "7" = "Inhib",
                              "8" = "OPCs",
                              "9" = "Ex",
                              "10" = "Oligos",
                              "11" = "Ex",
                              "12" = "Ex",
                              "13" = "Ex",
                              "14" = "Inhib",
                              "15" = "Ex",
                              "16" = "Inhib",
                              "17" = "Endo",
                              "18" = "Micro",
                              "19" = "Astros",
                              "20" = "Oligos",
                              "21" = "Inhib",
                              "22" = "Ex",
                              "23" = "Inhib",
                              "24" = "Ex",
                              "25" = "Astros",
                              "26" = "Ex",
                              "27" = "Oligos",
                              "28" = "Oligos")   ###new-brs

table(gse213.sc$cell.types)
color <- c(
  "Astros" = "#8DD3C7",
  "Ex" = "#FB8072",
  "Endo" = "#FFED6F",
  "Inhib" = "#80B1D3",
  "Micro" = "#BEBADA",
  "Oligos" = "#FDB462",
  "OPCs" = "#B3DE69")
p10 = DimPlot(gse213.sc, group.by = "cell.types", reduction = "umap", label = T, pt.size = 0.1, cols = color,raster=FALSE)
p10
ggsave(file.path(cellanno_dir, "09_CellType786452.pdf"), width = 10, height = 8, plot = p10)
p11 = DimPlot(gse213.sc, group.by = "cell.types", reduction = "umap", split.by = "group",
              ncol=2, label = F, pt.size = 0.1, cols = color,raster=FALSE)
p11
ggsave(file.path(cellanno_dir, "10_CellType786452_group.pdf"), width = 12, height = 5, plot = p11)

##
Idents(gse213.sc) <- gse213.sc$cell.types
markers_list <- list(
  Ex     = c("NRGN","TUBA1B","SLC17A7"),  #Excitatory Neurons(Ex) 
  Oligos = c("CERCAM","TF","ENPP2"),   # Oligodendrocytes(Oligos)
  Inhib  = c("GAD1","GAD2"), # Inhibitory Neurons (Inhib) 
  Astros = c("ADGRV1","SLC1A2","FGFR3"),  #Astrocytes (Astros) 
  Endo   = c("CLDN5","FLT1","ABCB1"),  # Endothelial cells (Endo)
  Micro  = c("P2RY12","CSF1R","DOCK8"),   # Microglia (Micro)
  OPCs   = c("PDGFRA","VCAN","SEMA5A")  # OPCs (OPCs) 
)
markers_ordered <- unlist(markers_list)
cell_type_order <- c(
  "Ex",
  "Oligos",
  "Inhib",
  "Astros",
  "Endo",
  "Micro",
  "OPCs"
)
table(gse213.sc$cell.types)
gse213.sc$cell_type_ordered <- factor(gse213.sc$cell.types, 
                                      levels = cell_type_order)

Idents(gse213.sc) <- "cell_type_ordered"
p9 <- DotPlot(
  object = gse213.sc,
  features = markers_list,                 
  group.by = "cell_type_ordered",      
  cols = c("lightgrey", "red"),       
  dot.scale = 6,                     
  assay = "RNA"
) +
  RotatedAxis() +                     
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),  
    plot.title.position = "plot"     
  ) +
  labs(title = "Marker Gene Expression by Cell Type")

print(p9)
ggsave(file.path(cellanno_dir, "10_CellType-marker.pdf"), width = 12, height = 6, plot = p9)
ggsave(file.path(cellanno_dir, "10_CellType-marker.tiff"), width = 12, height = 6, plot = p9, dpi = 300)

## 
library(dplyr)
library(tidyr)
Idents(gse213.sc) <- gse213.sc$orig.ident
gse213.sc$group <- as.factor(gse213.sc$group)
gse213.sc$Condition = ifelse(gse213.sc$group=='Case','MDD','CTL')
table(gse213.sc$sample, gse213.sc$group)
meta <- gse213.sc@meta.data %>%
  as.data.frame() %>%
  select(sample, cell.types, Condition)

Cell_type_stat_table <- meta %>%
  group_by(sample, cell.types, Condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = cell.types,
    values_from = count,
    values_fill = 0  # 缺失值填 0（如某个样本无 Microglia）
  ) %>%
  select(
    sample,
    Ex, Oligos, Astros, Inhib, OPCs, Endo, Micro,
    Condition
  )
print(Cell_type_stat_table)
write.csv(Cell_type_stat_table, file = file.path(cellanno_dir, "Cellnumber_count.csv"), row.names = FALSE)  
meta <- gse213.sc@meta.data %>%
  as.data.frame() %>%
  select(sample, group, Condition)
group_stat_table <- meta %>%
  group_by(sample, group, Condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = group,
    values_from = count,
    values_fill = 0  # 缺失值填 0
  ) %>%
  
  select(
    sample,
    Case,Control,
    Condition
  )
print(group_stat_table)
write.csv(group_stat_table, file = file.path(cellanno_dir, "sample.csv"), row.names = FALSE)   

#
library(ggthemes)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(ggprism)

table(gse213.sc$Condition)
cell_type_counts <- table(gse213.sc$cell.types, gse213.sc$Condition)
df <- as.data.frame(cell_type_counts) %>%
  setNames(c("cell_type", "condition", "count"))

proportion_df <- df %>%
  group_by(condition) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

condition_order <- c("CTL","MDD")
proportion_df$condition <- factor(proportion_df$condition, levels = condition_order, ordered = TRUE)
custom_colors <- c(
  "Astros" = "#8DD3C7",
  "Ex" = "#FB8072",
  "Endo" = "#FFED6F",
  "Inhib" = "#80B1D3",
  "Micro" = "#BEBADA",
  "Oligos" = "#FDB462",
  "OPCs" = "#B3DE69"
)
p11 <- ggplot(proportion_df, 
              aes(x = condition, 
                  y = percentage, 
                  fill = cell_type,
                  stratum = cell_type, 
                  alluvium = cell_type)) +
  
  geom_stratum(width = 0.7, color = 'white', alpha = 0.9) +
  geom_flow(alpha = 0.6, width = 0.7, color = 'white', size = 0.8, curve_type = "linear") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(mult = 0.3)) +
  labs(x = "Sample Type", y = "Relative Abundance (%)", fill = "Cell Type",
       title = "Cell Type Composition") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1)) +
  scale_fill_manual(values = custom_colors) +  
  theme(legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


plot(p11)
ggsave(file.path(cellanno_dir, "11_celltype_proportion_histogram.pdf"), width = 6, height = 8, plot = p11)
ggsave(file.path(cellanno_dir, "11_celltype_proportion_histogram.tiff"), width = 6, height = 8, dpi = 300, plot = p11)
save(gse213.sc, file = file.path(data_dir, "scRNA_brsrp.rdata"))

#########differential expression analysis######
set.seed(786452)
#######MAST#######
gse213.sc$class <- paste(gse213.sc$group, gse213.sc$cell.types, sep = "-")
Idents(gse213.sc) <- gse213.sc$class
library(MAST)
cell <- c('Astros','Endo','Ex','Inhib','Micro','Oligos','OPCs')
rt = gse213.sc@meta.data
a = unique(rt$class)
a = a[order(a)]
sc.marker = list()
for (i in 1:7) {
  marker1 <- FindMarkers(
    object = gse213.sc,
    ident.1 = a[i],           
    ident.2 = a[i + 7],       
    only.pos = FALSE,         
    logfc.threshold = 0.1,
    min.pct = 0.01,
    test.use = "MAST"
  )
  print(marker1[c("GLRA1","EFCAB2","ZNF605","TTC12","EVI5","DNAJC9","MZT2A"), , drop = FALSE])
  sc.marker[[i]] <- marker1
  print(i)
}
dataall <- list()
for (i in 1:7) {
  data <- sc.marker[[i]]
  data_df <- data.frame(
    gene_name = rownames(data),
    data,
    stringsAsFactors = FALSE
  )
  data_df$cell <- cell[i]
  dataall[[i]] <- data_df
}
dataall <- do.call(rbind, dataall)
rownames(dataall) <- NULL
write.csv(dataall,file.path(DEanalysis_dir, "allcell_DEGs.csv"))
xx <- c('Astros','Endo','Ex','Inhib',
        'Micro','Oligos','OPCs')
for (i in 1:7) {
  sc.marker1 = sc.marker[[i]]
  sc.marker1$cell = xx[i]
  write.csv(sc.marker1,file = paste0(file.path(DEanalysis_dir),xx[i],'_de.csv'))
}

#####single cell vocanol plot
library(ggplot2)
library(dplyr)
library(Matrix)
library(scales)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
multiVolcanoPlot = function(dat, color.arr=NULL, onlyAnnotateUp=T,
                            log2Foldchang=0.58, adjp=0.05, top_marker=5, 
                            max_overlaps=10, width=0.9){
  library(dplyr)
  library(ggrepel)
  if(is.null(color.arr)){
    len = length(unique(dat$cluster))
    color.arr=scales::hue_pal()(len)
  }
  
  dat.plot <- dat %>% mutate(
    "significance"=case_when(p_val_adj < adjp & avg_log2FC >= log2Foldchang  ~ 'Up',
                             p_val_adj < adjp & avg_log2FC <= -log2Foldchang  ~ 'Down',
                             TRUE ~ 'None'))
  tbl = table(dat.plot$significance)
  print( tbl )
  background.dat <- data.frame(
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC>0) %>%
      summarise("y.localup"=max(avg_log2FC)),
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC<=0) %>%
      summarise("y.localdown"=min(avg_log2FC)),
    x.local=seq(1:length(unique(dat.plot$cluster)))
  ) %>% dplyr::select(-cluster.1)
  
  x.number <- background.dat %>% dplyr::select(cluster, x.local)
  dat.plot <- dat.plot%>% left_join(x.number,by = "cluster")
  dat.marked.up <- dat.plot %>% filter(significance=="Up") %>%
    group_by(cluster) %>% arrange(-avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked.down <- dat.plot %>% filter(significance=="Down") %>%
    group_by(cluster) %>% arrange(avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked <- dat.marked.up %>% bind_rows(dat.marked.down)
  dat.infor <- background.dat %>%
    mutate("y.infor"=rep(0,length(cluster)))
  ##plotting:
  vol.plot <- ggplot()+
    geom_col(background.dat,mapping=aes(x.local, y.localup),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    geom_col(background.dat,mapping=aes(x.local,y.localdown),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    geom_jitter(dat.plot, mapping=aes(x.local, avg_log2FC, #x= should be number, Not string or factor
                                      color=significance),
                size=0.8, width = 0.4, alpha= 1)+
    scale_color_manual(name="significance", 
                       breaks = c('Up', 'None', 'Down'),
                       values = c("#d56e5e","#cccccc", "#5390b5")) + #set color for: Down None   Up
    geom_tile(dat.infor, mapping=aes(x.local, y.infor), #x axis color box
              height = log2Foldchang*1.3,
              fill = color.arr[1:length(unique(dat.plot$cluster))],
              alpha = 0.5,
              width=width) +
    labs(x=NULL,y="log2 Fold change")+
    geom_text(dat.infor, mapping=aes(x.local,y.infor,label=cluster))+
    ggrepel::geom_label_repel(data=if(onlyAnnotateUp) dat.marked.up else dat.marked, #gene symbol, of up group default
                              mapping=aes(x=x.local, y=avg_log2FC, label=gene),
                              force = 2, #size=2,
                              max.overlaps = 10000,
                              label.size = 0, #no border
                              fill="#00000000", #box fill color
                              seed = 111,
                              min.segment.length = 0,
                              force_pull = 2,
                              box.padding = 0,
                              segment.linetype = 1,
                              hjust = 0.5)+
    annotate("text", x=1.5, y=max(background.dat$y.localup)+3,
             label=paste0("|log2FC|>=", log2Foldchang, " & FDR<", adjp))+
    theme_classic(base_size = 12)+
    
    theme(
      axis.title = element_text(size = 13, color = "black"),
      axis.text = element_text(size = 15, color = "black"),
      axis.line.y = element_line(color = "black", size = 0.8),
      axis.line.x = element_blank(), #no x axis line
      axis.ticks.x = element_blank(), #no x axis ticks
      axis.title.x = element_blank(), #
      axis.text.x = element_blank(),
      legend.spacing.x = unit(0.1,'cm'),
      legend.key.width = unit(0.5,'cm'),
      legend.key.height = unit(0.5,'cm'),
      legend.background = element_blank(),
      legend.box = "horizontal",
      legend.position = c(0.08, 0.85),legend.justification = c(1,0)
    )+
    guides( #color = guide_legend( override.aes = list(size=5) ), #legend circle size
      color=guide_legend( override.aes = list(size=5), title="Change")
    )
  vol.plot
}
alldegs <- read.csv(file.path(DEanalysis_dir, "allcell_DEGs.csv"))
colors <- c("#8DD3C7","#FFED6F","#FB8072","#80B1D3",
            "#BEBADA","#FDB462","#B3DE69")
alldegs <- alldegs %>%
  mutate(label = case_when(
    p_val_adj < 0.05 & avg_log2FC >= 0.2  ~ "Up",
    p_val_adj < 0.05 & avg_log2FC <= -0.2 ~ "Down",
    TRUE ~ "None"
  ))
alldegs$label = factor(alldegs$label, levels = c("Up", "None", "Down"))
table(alldegs$label)
alldegs$cluster = factor(alldegs$cell, 
                         levels = c("Astros","Endo","Ex",
                                    "Inhib","Micro","Oligos","OPCs"))
alldegs <- alldegs %>%
  dplyr::rename(gene = gene_name)
newv <- multiVolcanoPlot(alldegs[,2:10],  
                         colors, onlyAnnotateUp = FALSE, top_marker = 3, 
                         log2Foldchang=0.2, adjp=0.05)
newv
ggsave(newv, file=file.path(DEanalysis_dir, "AllDEGs_valcano0.2_v1.pdf"),width =14, height = 6) 

####edgeR pseudobulk
rm(list = ls())
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead","#ffad60",
          "#7F95D1","#d9534f",  "#57D1C9", "#ED5485",  "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(BiocParallel))
library(patchwork)
####
load(file.path(data_dir, "scRNA_brsrp.rdata"))
Idents(gse213.sc) <- gse213.sc$orig.ident
gse213.sc$group <- as.factor(gse213.sc$group)
gse213.sc$Condition = ifelse(gse213.sc$group=='Case','MDD','CTL')
unique(gse213.sc$cell.types) ##查看所有细胞类型
Micro <- subset(gse213.sc,subset = cell.types == "Micro") ##根据细胞类型修改
###If Ex
###Ex <- subset(gse213.sc,subset = cell.types == "Ex")
######
run_pseudobulk_edger <- function(
    seu,
    sample_col = "sample",
    group_col = "Condition",
    min_cells = 10,
    output_file = NULL,
    robust = TRUE
){
  
  ## donor filter
  sample_vec <- seu@meta.data[[sample_col]]
  cellnum <- table(sample_vec)
  keep.sample <- names(
    cellnum[cellnum >= min_cells])
  cells.keep <- rownames(seu@meta.data
  )[sample_vec %in% keep.sample]
  seu <- subset(seu,cells = cells.keep)
  cat("Retained donors:",
      length(unique(seu[[sample_col]][,1])),"\n")
  
  ## pseudobulk
  counts <- GetAssayData(seu,assay = "RNA",layer = "counts")
  pb_counts <- rowsum(t(as.matrix(counts)),
                      group = seu[[sample_col]][,1])
  pb_counts <- t(pb_counts)
  
  ## sample metadata
  sample_info <- seu@meta.data %>%
    dplyr::select(all_of(sample_col),all_of(group_col)) %>%
    distinct()
  colnames(sample_info) <- c("sample", "Condition")
  
  ## sex info
  sample_info$sex <- ifelse(substr(sample_info$sample,1,1)=="F",
                            "Female","Male")
  sample_info$sex <- factor(sample_info$sex)
  
  ## 对齐
  sample_info <- sample_info[
    match(colnames(pb_counts),sample_info$sample),
  ]
  rownames(sample_info) <- sample_info$sample
  stopifnot(
    all(
      sample_info$sample ==
        colnames(pb_counts)
    )
  )
  
  ## edgeR
  dge <- DGEList(
    counts = pb_counts)
  keep <- filterByExpr(dge,
                       group = sample_info$Condition)
  dge <- dge[keep,,keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  sample_info$Condition <- factor(sample_info$Condition,
                                  levels = c("CTL","MDD"))
  design <- model.matrix(~ sex + Condition,
                         data = sample_info) ###添加性别作为协变量
  dge <- estimateDisp(dge,design)
  fit <- glmQLFit(dge,design,robust = robust)
  res <- glmQLFTest(fit,coef = "ConditionMDD")
  deg <- topTags(res,n = Inf)$table
  
  ## 保存
  if(!is.null(output_file)){
    write.csv(
      deg,
      output_file
    )
  }
  
  return(
    list(
      DEG = deg,
      fit = fit,
      dge = dge,
      design = design,
      sample_info = sample_info
    )
  )
}
micro_res <- run_pseudobulk_edger(seu = Micro,
                                  output_file =file.path(edgeR_dir, "Micro_edgeR.csv"))
##change cell Rdata name
#ex_res <- run_pseudobulk_edger(seu = Ex,
#                                  output_file =file.path(edgeR_dir, "Ex_edgeR.csv"))

####MAST-RE
rm(list = ls())
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead","#ffad60",
          "#7F95D1","#d9534f",  "#57D1C9", "#ED5485",  "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(BiocParallel))
suppressMessages(library(MAST))
library(patchwork)
####封装函数版
load(file.path(data_dir, "scRNA_brsrp.rdata"))
Idents(gse213.sc) <- gse213.sc$orig.ident
gse213.sc$group <- as.factor(gse213.sc$group)
gse213.sc$Condition = ifelse(gse213.sc$group=='Case','MDD','CTL')
unique(gse213.sc$cell.types) ##查看所有细胞类型
Micro <- subset(gse213.sc,subset = cell.types == "Micro") ##change based on celltypes
####RNA
# define and fit the model
run_MAST_RE_RNA <- function(
    seu,
    ref_group = "CTL",
    min_cells_per_sample = 10,
    min_pct = 0.1,
    out_file = NULL
){
  
  ## donor filter
  sample_n <- table(seu$sample)
  keep_samples <- names(
    sample_n[sample_n >= min_cells_per_sample])
  seu <- subset(seu, subset = sample %in% keep_samples)
  cat("Retained samples:", length(unique(seu$sample)), "\n")
  cat("Remaining cells:", ncol(seu), "\n\n")
  
  ## Add sex information
  if(!"sex" %in% colnames(seu@meta.data)){
    seu$sex <- ifelse(substr(seu$sample, 1, 1) == "F",
                      "Female", "Male")}
  seu$sex <- factor(seu$sex,
                    levels = c("Female","Male"))
  
  ##使用RNA assay
  DefaultAssay(seu) <- "RNA"
  if(ncol(GetAssayData(seu, layer = "data")) == 0){
    seu <- NormalizeData(seu,
                         normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE
    )
  }
  
  sce <- as.SingleCellExperiment(seu, assay = "RNA")
  sca <- SceToSingleCellAssay(sce)
  cat("Before gene filtering:\n")
  print(dim(sca))
  
  sca <- sca[freq(sca) > min_pct, ]
  cat("After gene filtering:\n")
  print(dim(sca))
  
  ##Cellular detection rate
  cdr2 <- colSums(assay(sca) > 0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  colData(sca)$Condition <- factor(
    colData(sca)$Condition)
  
  ## CTL as reference
  colData(sca)$Condition <- relevel(
    colData(sca)$Condition,
    ref = ref_group)
  colData(sca)$sex <- factor(
    colData(sca)$sex)
  colData(sca)$sample <- factor(
    colData(sca)$sample)
  
  ##MAST mixed model
  cat("\nFitting MAST model...\n")
  zlmCond <- zlm(formula = ~ Condition + sex + cngeneson + (1 | sample),
                 sca = sca, method = "glmer", ebayes = FALSE,
                 strictConvergence = FALSE, fitArgsD = list(nAGQ = 0)) 
  
  ##MDD vs CTL
  summaryCond <- summary(zlmCond, doLRT = "ConditionMDD")
  lrt.dt <- summaryCond$datatable
  result <- merge(lrt.dt[contrast=='ConditionMDD' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                  lrt.dt[contrast=='ConditionMDD' & component=='logFC', .(primerid, coef)],
                  by='primerid') # logFC coefficients
  ##lnFC -> log2FC
  result[, log2FC := coef / log(2)]
  ##FDR
  result[, FDR :=
           p.adjust(`Pr(>Chisq)`, method = "BH")]
  result <- as.data.frame(result)
  result <- result[
    order(result$FDR),]
  rownames(result) <- NULL
  ## save results
  write.csv(result, file = out_file, row.names = FALSE)
  cat( "\nResults saved to:", out_file, "\n")
  
  return(result)
}

mastRE_res <- run_MAST_RE_RNA(seu = Micro, #change based on celltypes
                              out_file = file.path(MASTRE_dir, "Micro_MAST_RE_RNA.csv"))

#####single-gene GSEA
library(Seurat)
library(corrplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(psych)
library(RColorBrewer)
library(GseaVis)
library(ggsci)
library(scales)
library(msigdbr)
library(limma)
library(psych)
library(DOSE)
load(file.path(data_dir, "scRNA_brsrp.rdata"))
set.seed(786452)
DimPlot(gse213.sc, group.by = "cell.types", reduction = "umap", label = T, pt.size = 0.1,raster=FALSE)
mycolor = pal_npg( "nrc", alpha = 0.8)(10) 
gse213.mdd = gse213.sc[,gse213.sc$Condition == "MDD"] 
table(gse213.mdd$cell.types)
gse213.mdd = gse213.mdd[,gse213.mdd$cell.types == "Ex",]  
Idents(gse213.mdd) = gse213.mdd$sample
expr = AverageExpression(gse213.mdd, assays = "SCT", slot = "data")[[1]]
dim(gse213.mdd)
dim(expr)
expr = expr[rowSums(expr)>0,]  
expr = data.frame(expr)
rt = expr
#
set.seed(3456)
output_dir = sgGSEA_dir
# Ex
gene = c("COQ8A")
options(timeout = 99999) 
for (ii in gene) {
  tryCatch({
    tar.exp = rt[ii, ]
    y = as.numeric(tar.exp)
    
    data1 = data.frame()
    for (i in rownames(rt)) {
      dd = corr.test(as.numeric(rt[i, ]), y, method = "spearman", adjust = "fdr")
      data1 = rbind(data1, data.frame(
        gene = i,
        cor = dd$r,
        p.value = dd$p,
        stringsAsFactors = FALSE
      ))
    }
    data1 = data1[order(data1$cor, decreasing = TRUE), ]
    
    gene1 = data1$cor
    names(gene1) = mapIds(
      org.Hs.eg.db,
      keys = data1$gene,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"  
    )
    table(is.na(names(gene1)))
    gene1 = gene1[!is.na(names(gene1))]
    
    if (length(gene1) == 0) {
      warning("无有效基因映射到 Entrez ID：", ii)
      next
    }
    
    #  GSEA
    gsea <- gseGO(
      geneList = gene1,
      OrgDb = org.Hs.eg.db,
      ont = "BP",keyType = "ENTREZID",
      pvalueCutoff = 1,verbose = T,seed = T)
    
    csv_path = file.path(output_dir, paste0(ii, '.Ex.3456.GOBP.csv'))  ###
    write.csv(as.data.frame(gsea), file = csv_path, row.names = FALSE)
    message("已保存：", csv_path)
    
    if (length(gsea) == 0 || nrow(as.data.frame(gsea)) == 0) {
      message("GSEA 无结果：", ii)
      next
    }
    
    gsea_df = as.data.frame(gsea)
    sig_ids = rownames(gsea_df[gsea_df$p.adjust < 0.05, ])  # FDR < 0.05
    sig_gsea = gsea_df[sig_ids,]
    csv_path = file.path(output_dir, paste0(ii, '.Ex.3456.GOBP.FDR0.05.csv')) ###
    write.csv(sig_gsea, file = csv_path, row.names = FALSE)
    message("已保存：", csv_path)
    if (length(sig_ids) == 0) {
      message("无显著富集通路 (FDR < 0.05)：", ii)
    } else {
      pdf_path = file.path(output_dir, paste0(ii, '.Ex.3456.GOBP.pdf'))   ###
      p1 = gseaplot2(
        gsea,                    
        geneSetID = sig_ids[1:min(10, length(sig_ids))],  
        subplots = 1:2,
        base_size = 16,
        title = paste0(ii, '-Ex'),   ####
        color = pal_d3(palette = "category20")(min(10, length(sig_ids)))
      )
      ggsave(pdf_path, width = 13, height = 12, plot = p1)
      message("已保存：", pdf_path)
    }
  }, error = function(e) {
    message("错误 [", ii, "]: ", e$message)
  })
}

#######
######COQ8A+ Ex analysis
rm(list = ls())
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead","#ffad60",
          "#7F95D1","#d9534f",  "#57D1C9", "#ED5485",  "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(harmony))
suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(BiocParallel))
library(patchwork)
load(file.path(data_dir, "scRNA_brsrp.rdata"))
ex <- subset(gse213.sc, subset = cell.types == "Ex")
DefaultAssay(ex) <- "SCT"
coq <- FetchData(ex, vars = "COQ8A", layer = "data")[,1]
ex$COQ8A_group <- ifelse(coq > 0, "COQ8A_pos", "COQ8A_neg")
table(ex$COQ8A_group)
###UMAP
gse213.sc$highlight <- "Other"
gse213.sc$highlight[Cells(ex)[ex$COQ8A_group == "COQ8A_pos"]] <- "COQ8A+ Ex"
set.seed(786452)
Idents(gse213.sc) = "SCT_snn_res.0.6"
gse213.sc$highlight <- factor(
  gse213.sc$highlight,
  levels = c("Other", "COQ8A+ Ex"))
####
gse213.sc$COQ8A_plot <- "Other cells"
idx_ex <- Cells(ex)
gse213.sc$COQ8A_plot[idx_ex] <- ex$COQ8A_group
gse213.sc$COQ8A_plot <- factor(
  gse213.sc$COQ8A_plot,
  levels = c("Other cells", "COQ8A_neg", "COQ8A_pos"))
p_umap_dis <- DimPlot(
  gse213.sc,
  reduction = "umap",        
  split.by = "Condition",
  group.by = "COQ8A_plot",
  cols = c("grey90", "grey80", "#D73027"),
  pt.size = 0.1,
  raster = FALSE)
p_umap_dis
ggsave(file.path(COQ8A_dir, "CellType_COQ8A_Ex_groupsplit.pdf"), width = 14, height = 8, p_umap_dis)
######
df_prop <- ex@meta.data %>%
  group_by(Condition) %>%
  summarise(
    total = n(),
    COQ8A_pos = sum(COQ8A_group == "COQ8A_pos"),
    prop = COQ8A_pos / total)
p2<-ggplot(df_prop, aes(x = Condition, y = prop, fill = Condition)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c(
      "MDD" = "#ea7070",
      "CTL" = "#005995")) +
  theme_classic() +
  labs(
    y = "Proportion of COQ8A+ Ex",
    x = NULL)
ggsave(file.path(COQ8A_dir, "Cellratio_COQ8A_Ex.pdf"), width = 6, height = 5, p2)
#########group set
ex$group_2x2 <- paste(ex$Condition, ex$COQ8A_group, sep = "_")
table(ex$group_2x2)
##########
#######
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(sctransform))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
###
deg_ctrl_ths <- FindMarkers(ex,
                            ident.1 = "CTL_COQ8A_pos",
                            ident.2 = "CTL_COQ8A_neg",
                            group.by = "group_2x2",
                            logfc.threshold = 0.26,
                            min.pct = 0.05)
deg_dis_ths <- FindMarkers(ex,
                           ident.1 = "MDD_COQ8A_pos",
                           ident.2 = "MDD_COQ8A_neg",
                           group.by = "group_2x2",
                           logfc.threshold = 0.26,
                           min.pct = 0.05)
write.csv(deg_ctrl_ths, file = file.path(COQ8A_dir, "COQ8A_CTL_DEGs_0.26.csv"))
write.csv(deg_dis_ths, file = file.path(COQ8A_dir, "COQ8A_MDD_DEGs_0.26.csv"))
###
deg_ctrl <- FindMarkers(ex,
                        ident.1 = "CTL_COQ8A_pos",
                        ident.2 = "CTL_COQ8A_neg",
                        group.by = "group_2x2",
                        logfc.threshold = 0,
                        min.pct = 0.05)
deg_dis <- FindMarkers(ex,
                       ident.1 = "MDD_COQ8A_pos",
                       ident.2 = "MDD_COQ8A_neg",
                       group.by = "group_2x2",
                       logfc.threshold = 0,
                       min.pct = 0.05)
write.csv(deg_ctrl, file = file.path(COQ8A_dir, "COQ8A_CTL_DEGs.csv"))
write.csv(deg_dis, file = file.path(COQ8A_dir, "COQ8A_MDD_DEGs.csv"))

###
#######COQ8A+ Ex vocanol plots
library(tools)
library(dplyr)
library(readr)
library(ggvenn)
library(ggrastr)
########
#res_8AMDD <- read.csv(file.path(COQ8A_dir, "COQ8A_MDD_DEGs.csv"))
res_8AMDD <- read.csv(file.path(COQ8A_dir, "COQ8A_CTL_DEGs.csv")) ##res_8ACTL
#######
res_8AMDD$type <- case_when(res_8AMDD$avg_log2FC > 0.26 & res_8AMDD$p_val_adj < 0.05 ~ 'Up',
                            res_8AMDD$avg_log2FC < -0.26 & res_8AMDD$p_val_adj < 0.05 ~ 'Down',
                            TRUE ~ 'None')

res_8AMDD$p_val_adj <- as.numeric(as.character(res_8AMDD$p_val_adj))
res_8AMDD$p_val_adj[res_8AMDD$p_val_adj == 0] <- 1e-300  ###赋值
res_8AMDD <- res_8AMDD %>% 
  filter(!is.na(avg_log2FC), !is.na(p_val_adj))
res_8AMDD1 <- res_8AMDD %>% filter(avg_log2FC < 10 & avg_log2FC > -10)
cut_off_padj = 0.05
cut_off_logFC = 0.26
p1 <- ggplot(res_8AMDD1, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=type)) +
  rasterise(geom_point(alpha=0.8, size=1), dpi = 300) +
  scale_color_manual(values=c("#005995", "#d2dae2","#d73027"))+
  # 辅助线
  geom_vline(xintercept=c(-0.26,0.26),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_padj),
             lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(FoldChange)",
       y="-log10(padjust)")+
  theme_bw()+
  ggtitle("Volcano Plot of CTL COQ8A+/COQ8A- Ex DEGs")+  #MDD
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

p2 <- p1 + coord_flip()
p2
ggsave(file.path(COQ8A_dir, "CTL_COQ8A_volcano_simple_withoutEX.pdf"), p2, device = cairo_pdf, width = 7, height = 6)

#####GO
library(tools)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
folder_path <- COQ8A_dir
file_suffix <- "_DEGs_0.26.csv"
file_list <- list.files(
  path = folder_path,
  pattern = paste0(".*", file_suffix, "$"), 
  full.names = TRUE,                         
  ignore.case = FALSE                        
)
file_list
fixed_prefix <- COQ8A_dir
fixed_suffix <- "_DEGs_0.26.csv"
options(timeout = 99999) 
for (file in file_list) {
  cat("\n==== 正在处理文件:", basename(file), "====\n")
  data <- read.csv(file)
  data <- data[data$X != "COQ8A", ] 
  if ("p_val_adj" %in% colnames(data) && "avg_log2FC" %in% colnames(data)) {
    sig_genes <- data %>% 
      filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.26)} 
  gene_symbols <- sig_genes[, 1]
  genelist <- bitr(gene_symbols, fromType="SYMBOL",
                   toType="ENTREZID", OrgDb='org.Hs.eg.db')
  genelist <- pull(genelist,ENTREZID)
  ##GO
  ego_ALL <- enrichGO(gene = genelist, OrgDb = 'org.Hs.eg.db',
                      ont = "ALL",  
                      pAdjustMethod = "BH", pvalueCutoff  = 0.05,readable = TRUE)
  middle_part <- sub(paste0("^", fixed_prefix, "(.*)", fixed_suffix, "$"), "\\1", file)
  write.csv(ego_ALL@result, file.path(folder_path, paste0(middle_part, "_enrichGO_padj.csv")))
}

############################GSEA
library(corrplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(psych)
library(GseaVis)
library(msigdbr)
library(limma)
library(genekitr)
library(data.table)
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
output_dir <- COQ8A_dir
deg_dis <- read.csv(file.path(COQ8A_dir, "COQ8A_MDD_DEGs.csv"))
deg_ctrl <- read.csv(file.path(COQ8A_dir, "COQ8A_CTL_DEGs.csv"))
set.seed(123456789)
options(timeout = 99999) 
## SYMBOL + FC
deg_ctrl_filt <- deg_ctrl[rownames(deg_ctrl) != "COQ8A", ]
deg_ctrl_filt$gene <- row.names(deg_ctrl_filt)
d <- deg_ctrl_filt %>%
  dplyr::select(SYMBOL = gene, FC = avg_log2FC) %>%
  dplyr::filter(!is.na(SYMBOL), FC != 0)
entrez <- bitr(d$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
final <- inner_join(d, entrez, by = "SYMBOL")
final2 <- final %>%
  group_by(ENTREZID) %>%
  summarise(FC = max(FC, na.rm = TRUE)) %>%
  arrange(desc(FC))
if (nrow(final2) < 10) {
  message("Too few ENTREZ genes for: ", i)
  return(NULL)}
genelist <- final2$FC
names(genelist) <- final2$ENTREZID
## ===== GO BP =====
gse_go <- gseGO(
  geneList = genelist,
  OrgDb = org.Hs.eg.db,
  ont = "BP",keyType = "ENTREZID",
  pvalueCutoff = 1,verbose = FALSE)
if (!is.null(gse_go) && nrow(gse_go@result) > 0) {
  fwrite(as.data.frame(gse_go@result),
         file = file.path(output_dir, paste0("COQ8A_CTL_GOBP_GSEA.csv")))}

####enrichment visualization
####GOBP & GSEA
#GOBP
library(dplyr)
library(tidyr)
library(ggplot2)
go_result <- read.csv(file.path(COQ8A_dir, "COQ8A_CTL_enrichGO_padj.csv"), row.names = 1, stringsAsFactors = FALSE)
#######
top30_go <- go_result %>%
  filter(ONTOLOGY == "BP") %>%
  mutate(logpadj = -log10(p.adjust)) %>%
  arrange(desc(logpadj)) %>%
  slice(1:30)
top30_go$Description <- factor(top30_go$Description, levels = rev(top30_go$Description))
######
max_logpadj <- max(top30_go$logpadj) 
max_count <- max(top30_go$Count) 
scaled_max_y <- max_logpadj * 1.1 
##########
p5 <- ggplot() +
  geom_bar(data = top30_go,
           aes(x = Description, y = logpadj, fill = logpadj),
           stat = "identity", width = 0.8) +
  
  scale_y_continuous(
    name = "-log10(p.adjust)",
    limits = c(0, scaled_max_y),  # 注意放大 y 的上限
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * max_count / max_logpadj,
                        name = "Count",
                        breaks = pretty(c(0, max_count), 5))
  )+
  
  geom_line(data = top30_go,
            aes(x = Description, y = Count * max_logpadj / max_count, group = 1),
            color = "grey", linewidth = 0.8) +
  geom_point(data = top30_go,
             aes(x = Description, y = Count * max_logpadj / max_count),
             size = 3) +
  
  scale_fill_gradientn(
    colors = c("#FDE0DD", "#FA9FB5", "#C51B8A", "#7A0177"))+
  labs(title = "GO-BP enrichment barplot", x = "", fill = "-log10(p.adjust)") +
  coord_flip(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.x.top = element_text(size = 14),
    legend.position = "right"
  )
p5
ggsave(file.path(COQ8A_dir, "CTL_COQ8AEx_GOtop30_bar_dotplot.pdf"), p5, width = 12, height = 8)

##GSEA
gsea_df <- read.csv(file.path(COQ8A_dir, "COQ8A_MDD_GOBP_GSEA.csv"))
library(dplyr)
gsea_df <- gsea_df %>% 
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)
gsea_df <- as.data.frame(gsea_df) %>%
  mutate(
    Count = sapply(strsplit(core_enrichment, "/"), length))
gsea_plot_df <- gsea_df %>%
  mutate(group = ifelse(NES > 0, "activated", "suppressed")) %>%
  group_by(group) %>%
  arrange(group,
          ifelse(group == "activated", -NES, NES)) %>%
  slice_head(n = 15) %>%
  ungroup()
gsea_plot_df <- gsea_plot_df %>%
  mutate(
    Description = factor(
      Description,
      levels = Description[order(NES)]
    )
  )
gsea_plot_df <- gsea_plot_df %>%
  arrange(NES) %>%   
  mutate(
    Description = factor(Description, levels = unique(Description))
  )
library(ggplot2)
p <- ggplot(data = gsea_plot_df,aes(x = NES, y = Description)) +
  geom_segment(aes(x = 0,xend = NES,y = Description,yend = Description))+
  geom_point(aes(color = p.adjust, size = Count))+
  scale_color_gsea()+
  scale_y_discrete(labels=function(x) str_wrap(x, width=32)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12, colour = "gray30")) +
  labs(y = NULL)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=35)) + 
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  labs(colour = "p.adjust")+labs(size = "Count")
p
p_final <- p +
  scale_color_gradientn(
    colors = c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9","#2C7BB6" )
  )

p_final
ggsave(file.path(COQ8A_dir, "MDD_COQ8AEx_GOBP_GSEA_lollipop.pdf"), p_final, width = 10, height = 8)

####supplemental analysis
#####edgeR NB-GLM
########Differential Abundance analysis
library(dplyr)
library(edgeR)
meta <- gse213.sc@meta.data
count_mat <- table(meta$cell.types,meta$sample)
count_mat <- as.matrix(count_mat)
sample_info <- meta[, c("sample", "Condition"), drop = FALSE]
sample_info <- unique(sample_info)
rownames(sample_info) <- sample_info$sample
sample_info$sex <- ifelse(
  grepl("^F", rownames(sample_info)), "Female", "Male")
sample_info$sample <- NULL
sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c("CTL", "MDD"))
sample_info$sex <- factor(
  sample_info$sex,
  levels = c("Male", "Female")  
)
###check
table(sample_info$Condition, sample_info$sex)
######
dge <- DGEList(
  counts = count_mat,
  samples = sample_info)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ Condition + sex, data = dge$samples) ##加入sex作为协变量
colnames(design)
####NB-GLM
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
res <- glmQLFTest(fit, coef = "ConditionMDD")
da_res <- topTags(res, n = Inf)$table
da_res$celltype <- rownames(da_res)
head(da_res)
write.csv(da_res, file = file.path(DAanalysis_dir, "DAresults_all.csv"),row.names = T)
sig_da <- da_res %>%
  filter(FDR < 0.05)
write.csv(sig_da, file = file.path(DAanalysis_dir, "DAresults_sig.csv"),row.names = T)
