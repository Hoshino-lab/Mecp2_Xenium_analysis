library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

xenium_paths <- c(
  "Region_1" = "output-XETG00310__0046278__Region_1__20250515__065552",
  "Region_2" = "output-XETG00310__0046278__Region_2__20250515__065552",
  "Region_7" = "output-XETG00310__0046278__Region_7__20250515__065552",
  "Region_8" = "output-XETG00310__0046278__Region_8__20250515__065552",
  "Region_5" = "output-XETG00310__0046292__Region_5__20250515__065552",
  "Region_6" = "output-XETG00310__0046292__Region_6__20250515__065552"
)
cell_id_paths <- c(
  "Ctrl_1" = "PFC_Control_1.csv",
  "Ctrl_2" = "PFC_Control_2.csv",
  "KO_1"   = "PFC_KO_1.csv",
  "KO_2"   = "PFC_KO_2.csv",
  "Rescue_1" = "PFC_Rescue_1.csv",
  "Rescue_2" = "PFC_Rescue_2.csv"
)


### Load Xenium output as Xenium object
region1_full_obj <- LoadXenium(data.dir = xenium_paths["Region_1"], fov = "Region_1_full")
region2_obj <- LoadXenium(data.dir = xenium_paths["Region_2"], fov = "Region_2_Ctrl_3") # WT_3
region7_full_obj <- LoadXenium(data.dir = xenium_paths["Region_7"], fov = "Region_7_full")
region8_obj <- LoadXenium(data.dir = xenium_paths["Region_8"], fov = "Region_8_KO_3")   # KO_3
region5_full_obj <- LoadXenium(data.dir = xenium_paths["Region_5"], fov = "Region_5_full")
region6_obj <- LoadXenium(data.dir = xenium_paths["Region_6"], fov = "Region_6_Rescue_3") # WT+3TC_3


### Cell ID

cell_ids_ctrl_1 <- read.csv(cell_id_paths["Ctrl_1"], header = TRUE)$Cell.ID
cell_ids_ctrl_2 <- read.csv(cell_id_paths["Ctrl_2"], header = TRUE)$Cell.ID
cell_ids_ko_1 <- read.csv(cell_id_paths["KO_1"], header = TRUE)$Cell.ID
cell_ids_ko_2 <- read.csv(cell_id_paths["KO_2"], header = TRUE)$Cell.ID
cell_ids_rescue_1 <- read.csv(cell_id_paths["Rescue_1"], header = TRUE)$Cell.ID
cell_ids_rescue_2 <- read.csv(cell_id_paths["Rescue_2"], header = TRUE)$Cell.ID


### Splitting object into individual animals
# --- Control (N=3) ---
# Region_1 -> WT-1, WT-2
ctrl_1_obj <- subset(region1_full_obj, cells = cell_ids_ctrl_1)
ctrl_2_obj <- subset(region1_full_obj, cells = cell_ids_ctrl_2)
rm(region1_full_obj) 

# Region_2 -> WT-3
ctrl_3_obj <- region2_obj

# --- KO (N=3) ---
# Region_7 ->KO_1 , KO_2 
ko_1_obj <- subset(region7_full_obj, cells = cell_ids_ko_1)
ko_2_obj <- subset(region7_full_obj, cells = cell_ids_ko_2)

# Region_8 -> KO_3
ko_3_obj <- region8_obj

# --- KO+3TC (N=3) ---
# Region_5 -> KO+3TC-1, KO+3TC-2 
rescue_1_obj <- subset(region5_full_obj, cells = cell_ids_rescue_1)
rescue_2_obj <- subset(region5_full_obj, cells = cell_ids_rescue_2)


# Region_6 -> KO+3TC-3
rescue_3_obj <- region6_obj
rm(region6_obj)

#adding metadata
# Control
ctrl_1_obj$individual <- "Ctrl_1"
ctrl_1_obj$condition <- "Control"
ctrl_2_obj$individual <- "Ctrl_2"
ctrl_2_obj$condition <- "Control"
ctrl_3_obj$individual <- "Ctrl_3"
ctrl_3_obj$condition <- "Control"

# KO
ko_1_obj$individual <- "KO_1"
ko_1_obj$condition <- "KO"
ko_2_obj$individual <- "KO_2"
ko_2_obj$condition <- "KO"
ko_3_obj$individual <- "KO_3"
ko_3_obj$condition <- "KO"

# Rescue
rescue_1_obj$individual <- "Rescue_1"
rescue_1_obj$condition <- "Rescue"
rescue_2_obj$individual <- "Rescue_2"
rescue_2_obj$condition <- "Rescue"
rescue_3_obj$individual <- "Rescue_3"
rescue_3_obj$condition <- "Rescue"


#integrating objects

all_individual_objects_list <- list(
  ctrl_1_obj, ctrl_2_obj, ctrl_3_obj,
  ko_1_obj, ko_2_obj, ko_3_obj,
  rescue_1_obj, rescue_2_obj, rescue_3_obj
)

merged.obj <- all_individual_objects_list[[1]]
for (i in 2:length(all_individual_objects_list)) {
  cat("Merging object", i, "of", length(all_individual_objects_list), "\n")
  merged.obj <- merge(merged.obj, y = all_individual_objects_list[[i]])
}


##Subsetting lowquality cells, Normalization, Dimensionality reduction

merged.obj <- subset(merged.obj, subset = nFeature_Xenium > 250 & nCount_Xenium > 350)
merged.obj <- NormalizeData(merged.obj)
merged.obj <- FindVariableFeatures(merged.obj)
merged.obj <- ScaleData(merged.obj, features = rownames(merged.obj))
merged.obj <- RunPCA(merged.obj, features = VariableFeatures(object = merged.obj))
ElbowPlot(merged.obj) 
dims_to_use <- 1:20
merged.obj <- RunUMAP(merged.obj, dims = dims_to_use)
merged.obj <- FindNeighbors(merged.obj, dims = dims_to_use)
merged.obj <- FindClusters(merged.obj, resolution = 1.0)

##Annotation
cell_type <- c("Excitatory neuron_1","Excitatory neuron_2","Excitatory neuron_3","Endothelial cell_1","Excitatory neuron_4",
               "Inhibitory neuron_1","Astrocyte_1","Inhibitory neuron_2","Excitatory neuron_5","Excitatory neuron_6",
               "Oligodendrocyte_1","Microglia_1","Excitatory neuron_7","Astrocyte_2","Excitatory neuron_8","Oligodendrocyte_2",
               "OPC_1","Fibroblast","Inhibitory neuron_3","Excitatory neuron_9","Smooth muscle cell","Excitatory neuron_10",
               "Inhibitory neuron_4","Inhibitory neuron_5","Excitatory neuron_11","Inhibitory neuron_6","Pericyte","Excitatory neuron_12",
               "Excitatory neuron_13","Oligodendrocyte_3","Oligodendrocyte_4","Astrocyte_3","Endothelial cell_2","Oligodendrocyte_4","Microglia_2","OPC_2")
names(cell_type) <- levels(merged.obj)
merged.obj <- RenameIdents(merged.obj, cell_type)
merged.obj$cell_type_detail <- Idents(merged.obj)
merged.obj$cell_type <- str_split(merged.obj$cell_type_detail,pattern = "_",simplify = TRUE)[,1]
Idents(merged.obj) <- "cell_type"

saveRDS(merged.obj, "wholePFC_named_final.rds")

Ex <- subset(so,idents = "Excitatory neuron")
Ex %<>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5) 

AS <- subset(so,idents = "Astrocyte")
AS %<>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5) 

saveRDS(AS, "Astrocyte_only.rds")
saveRDS(Ex, "Exicitatory_only.rds")