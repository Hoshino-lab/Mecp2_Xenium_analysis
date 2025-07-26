library(Seurat)
library(dplyr)

#gene list to compare

genes_to_compare <- c(
  "Abl1","Adar","Aldoa","Atp5f1","Brd2","Cbx5","Cltc","Creb1","Csde1",
  "Cyfip2","Dpysl2","Flna","Hnrnpdl","Hnrnph1","Hspa1b","Mcl1","Msn",
  "Myh9","Pkm","Polr2a","Prkar1a","Prrc2b","Rbbp4","Set",
  "Setd1b","Smarca4","Sptan1","Txnip","Ywhaz","Zeb2"
)

# Loading Astrocyte and Excitatory neuron

astrocyte.obj <- readRDS("Astrocyte_only.rds")
excitatory_neuron.obj <- readRDS("Excitatory_only.rds")
genes_common <- intersect(genes_to_compare, rownames(astrocyte.obj))

##Astrocyte Fisher's test

#Creating astrocyte DEG list
Idents(astrocyte.obj) <- "condition"
all_markers_astro <- FindAllMarkers(astrocyte.obj)
up_in_ko_astro <- subset(all_markers_astro, cluster == "KO" & p_val_adj < 0.05 & avg_log2FC > 0)
total_up_degs_astro <- nrow(up_in_ko_astro)

list_up_degs_astro <- subset(up_in_ko_astro, gene %in% genes_exist_astro)
num_list_up_degs_astro <- nrow(list_up_degs_astro)

#Making contingency table and test
total_genes_astro <- nrow(astrocyte.obj)
a <- num_list_up_degs_astro; b <- total_up_degs_astro - a; c <- length(genes_exist_astro) - a; d <- total_genes_astro - total_up_degs_astro - c
table_up_astro <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, dimnames = list(c("Is KO-UP", "Not KO-UP"), c("In List", "Not In List")))
print(" (Astrocyte, UP in KO vs Others):")
print(table_up_astro)
fisher_up_astro <- fisher.test(table_up_astro, alternative = "greater")
print(fisher_up_astro)

#Creating excitatory DEG list
Idents(excitatory_neuron.obj) <- "condition"
all_markers_excit <- FindAllMarkers(excitatory_neuron.obj)

#Filtering out gene names commonly observed in 'DEGs' and 'genes_to_compare'
up_in_ko_excit <- subset(all_markers_excit, cluster == "KO" & p_val_adj < 0.05 & avg_log2FC > 0)
total_up_degs_excit <- nrow(up_in_ko_excit)

list_up_degs_excit <- subset(up_in_ko_excit, gene %in% genes_exist_excit)
num_list_up_degs_excit <- nrow(list_up_degs_excit)

#Making contingency table and test

total_genes_excit <- nrow(excitatory_neuron.obj)
a <- num_list_up_degs_excit; b <- total_up_degs_excit - a; c <- length(genes_exist_excit) - a; d <- total_genes_excit - total_up_degs_excit - c
table_up_excit <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, dimnames = list(c("Is KO-UP", "Not KO-UP"), c("In List", "Not In List")))
print(" (Excitatory Neuron, UP in KO vs Others):")
print(table_up_excit)
fisher_up_excit <- fisher.test(table_up_excit, alternative = "greater")
print(fisher_up_excit)